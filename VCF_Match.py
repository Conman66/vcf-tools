#!/usr/bin/python

############################################################
#
# VCF_Match.py
#
# Uses Python 2.6.6
#
# Searches an arbitrary number of VCF files and finds 
# the genomic positions where all files match in their
# called allele and are also different than the reference.
#
# Author: Connor Howington
# Date: 6/28/17
#
############################################################

#---------------------------------------------------------------------------------------
# This script uses a nested 3D dictionary to build a unified view of the variations
# in the submitted VCF files.  This data structure is essentially a set of trees,
# each one representing the variations across all VCFs for a single chromosome.
#
# Instantiation: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0, []])))
# This sets the default value of the innermost dictionary to [0, []].
#
# Diagram:
#
# 	Dict:  Key    Value
#              CHROM  Dict:  Key  Value
#                            POS  Dict:  Key        Value
#                                        "REF/ALT"  List
#                                                   [#, List]
#                                                       [QUALs]
#
# CHROME and POS are the same as in VCF, "REF/ALT" is a string containing REF and ALT
# separated by "/", # is the number of VCFs that contain this variation, and the list
# of QUALs contains the quality that each VCF reported for that variation.
#
# Note: This structure makes efficient use of memory and allows for compact code,
# but the code is somewhat hard to read.  It might be worthwhile to implement this
# in a different way, possibly using classes.
#---------------------------------------------------------------------------------------

import argparse, os, sys
from collections import defaultdict


def parse_args():
	parser = argparse.ArgumentParser(description='Searches an arbitrary number of VCF files and finds the genomic positions where the files match in their called allele and are also different than the reference.  Output is printed to stdout.')
	parser.add_argument('VCF_file', help='an input VCF file', nargs='+')
	parser.add_argument('--remove', '-r', help='remove any variations that are also present within the given file(s)', dest='remove_file', nargs='+')

	mutex_group = parser.add_mutually_exclusive_group()
	mutex_group.add_argument('--all', '-a', action='store_true', help='display all variations, not just the ones that are shared across all VCFs') 
	mutex_group.add_argument('--low_bound', '-l', type=int, help='only return variations that have occurred in at least LOW_BOUND VCF files')

	args = parser.parse_args()
	paths = list(args.VCF_file)
	
	if args.remove_file:
		paths += args.remove_file

	if args.low_bound and args.low_bound < 0:
		parser.error('Low bound must be at least 0')

	for path in paths:
		if not os.path.isfile(path):
			parser.error('File "{0}" cannot be found'.format(path))

	return args


# Takes a list of VCF file paths and returns a 3D dict representing a unified view
# of all the variations within those files
def find_variations(vcf_paths):
	variations = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0, []])))

	for path in vcf_paths:
		vcf_file = open(path, 'r')
		
		line = vcf_file.readline()

		# Skip meta-information
		while line[0:2] == '##':
			line = vcf_file.readline()

		# Check for header line
		if line[0] != '#':
			sys.stderr.write('Error in "{0}": Invalid or missing header line\n'.format(path))
			sys.exit(1)

		line = line.rstrip()
		split_line = line.split('\t')
		split_line[0] = split_line[0][1:]  # Remove '#' from first field

		# Find the indexes of all required fields based on header line
		# Exit with error if not all are found
		try:
			chrom_index = split_line.index('CHROM')
			pos_index = split_line.index('POS')
			ref_index = split_line.index('REF')
			alt_index = split_line.index('ALT')
			qual_index = split_line.index('QUAL')
		except ValueError:
			sys.stderr.write('Error in "{0}": Required fields missing in header line\n'.format(path))
			sys.exit(1)

		# Optional field
		#try:
		#	num_index = split_line.index('NUM')
		#except ValueError:
		#	num_index = None

		line = vcf_file.readline()

		# Iterate over all data lines
		while line != '':
			split_line = line.split('\t')
			
			# If a meta-info line is found within the data, exit with error
			if split_line[0] == '#':
				sys.stderr.write('Error in "{0}": Unexpected meta-information line within data\n'.format(path))
				sys.exit(1) 
			
			alt = split_line[alt_index]

			# If there is a variation
			if alt not in ['.', 'N']:
				chrom = split_line[chrom_index]
				pos = int(split_line[pos_index])
				ref = split_line[ref_index]
				qual = split_line[qual_index].rstrip().split(';')

				#if num_index:
				#	num = int(split_line[num_index])
				#else:
				num = 1

				ref_alt = ref + '/' + alt   # A string of the form "[ref]/[alt]" 
				
				# Increment the counter and add the reported quality to the quality list
				variations[chrom][pos][ref_alt][0] += num
				variations[chrom][pos][ref_alt][1] += qual

			line = vcf_file.readline()		

		vcf_file.close()

	return variations


# Takes a 3D dict containing all variations and returns a 3D dict containing only
# the variations that were shared across at least [low_bound] VCF files
def shared_variations(vars, low_bound):
	shared = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0, []])))

	for chrom in vars:
		for pos in vars[chrom]:
			for ref_alt in vars[chrom][pos]:
				if vars[chrom][pos][ref_alt][0] >= low_bound:
					shared[chrom][pos][ref_alt] = vars[chrom][pos][ref_alt]
	
	return shared


# Takes a 3D dict containing variations and prints each variation
# line-by-line in a VCF-like format
def print_variations(vars, num_files):
	print('#CHROM\tPOS\tREF\tALT\tNUM\tPER\tQUAL')

	for chrom in sorted(vars.keys()):
		for pos in sorted(vars[chrom].keys()):
			for ref_alt in vars[chrom][pos]:
				split_string = ref_alt.split('/')

				ref = split_string[0]
				alt = split_string[1]
				num = vars[chrom][pos][ref_alt][0]
				per = (float(num) / num_files) * 100
				quals = vars[chrom][pos][ref_alt][1]

				qual_string = quals[0]
				for i in range(1, len(quals)):
					qual_string += ';' + quals[i]

				print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}%\t{6}'.format(chrom, pos, ref, alt, num, round(per, 1), qual_string))

# Takes two 3D dicts and returns a 3D dict containing only the variations contained
# in the first dict and not the second (in other words, a set subtraction)
def remove_variations(source_vars, remove_vars):
	new_vars = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0, []])))

	for chrom in source_vars:
		for pos in source_vars[chrom]:
			for ref_alt in source_vars[chrom][pos]:
				if not ((chrom in remove_vars) and (pos in remove_vars[chrom]) and (ref_alt in remove_vars[chrom][pos])):
					new_vars[chrom][pos][ref_alt] = source_vars[chrom][pos][ref_alt]
	return new_vars


# Main flow
args = parse_args()
add_files = args.VCF_file
remove_files = args.remove_file
variations = find_variations(add_files)

if not args.all:
	if args.low_bound:
		low_bound = args.low_bound
	else:
		low_bound = len(add_files)

	variations = shared_variations(variations, low_bound)

if remove_files:
	remove_vars = find_variations(remove_files)
	variations = remove_variations(variations, remove_vars)

print_variations(variations, len(add_files))

