#!/usr/bin/python

#####################################################
#
# 1d_map_vcf.py
# 
# Uses a delta file created by MUMmer to map SNPs
# from one namespace to another in a VCF file
#
# Author: Connor Howington
# Date: 7/7/17
#
# Uses Python 2.6.6
#
#####################################################

import argparse, os, sys

def parse_args():
	parser = argparse.ArgumentParser(description='Uses a delta file created by MUMmer to map SNPs from one namespace to another in a VCF file ')

	parser.add_argument('VCF_file')
	parser.add_argument('delta_file')
	
	parser.add_argument('--reverse', '-r', help='Reverses the mapping process so that the reference space is translated into the query space', action='store_true')

	return parser.parse_args()


def create_mapping(delta_path, reverse):
	mapping = dict()
	gapmap = dict()

	delta_file = open(delta_path, 'r')
	line = delta_file.readline()

	while line[0] != '>' and line != '' :
		line = delta_file.readline()

	while line != '':
		#print(line)
		split_line = line.split()	

		if line[0] == '>':
			ref = split_line[0][1:]
			query = split_line[1]
			line = delta_file.readline()
			split_line = line.split()
		#print(line)
		ref_start = int(split_line[0])
		ref_end = int(split_line[1])
		query_start = int(split_line[2])
		query_end = int(split_line[3])

		#print(ref_start)

		# Determine directions so we can handle forward and backward strand matches
		ref_dir = cmp(ref_end - ref_start, 0)
		query_dir = cmp(query_end - query_start, 0)

		ref_pos = ref_start - ref_dir
		query_pos = query_start - query_dir

		num = delta_file.readline().rstrip()

		#abs_sum = 0
		
		last_match = None

		while num != '0':
			num = int(num)
			#abs_sum += abs(num)

			for i in range(1, abs(num)):
				ref_pos += ref_dir
				query_pos += query_dir
				
				if not reverse:
					mapping[query+'\t'+str(query_pos)] = ref+'\t'+str(ref_pos)
					last_match = ref_pos
				else:
					mapping[ref+'\t'+str(ref_pos)] = query+'\t'+str(query_pos)
					last_match = query_pos

			if num > 0:
				ref_pos += ref_dir
				if reverse and last_match is not None:
					gapmap[ref+'\t'+str(ref_pos)] = query+'\t'+str(last_match)
			else:
				query_pos += query_dir
				if not reverse and last_match is not None:
					gapmap[query+'\t'+str(query_pos)] = ref+'\t'+str(last_match)

			num = delta_file.readline().rstrip()
			#print(num)

		# If the current positions are not the same distance from the end, something's wrong with the file
		if (abs(ref_end - ref_pos)) != (abs(query_end - query_pos)):
			print('Error: Invalid delta file: Mismatched sequence numbers at\n{0} {1}\n{2} {3} {4} {5}'.format(ref, query, ref_start, ref_end, query_start, query_end))
			sys.exit(1)
		
		for i in range(abs(ref_end - ref_pos)):
			ref_pos += ref_dir
			query_pos += query_dir

			if not reverse:
				mapping[query+'\t'+str(query_pos)] = ref+'\t'+str(ref_pos)
			else:
				mapping[ref+'\t'+str(ref_pos)] = query+'\t'+str(query_pos)

		#print(abs_sum)
		#print(ref_pos)
		#print(query_pos)
		
		line = delta_file.readline()

	delta_file.close()
	return mapping, gapmap


def translate_vcf(vcf_path, mapping, gapmap):
	old_vcf = open(vcf_path, 'r')
	new_vcf = open(vcf_path+'.map', 'w', 1)
	gap_vcf = open(vcf_path+'.gap', 'w', 1)

	line = old_vcf.readline()

	while line[0] == '#':
		new_vcf.write(line)
		gap_vcf.write(line)
		line = old_vcf.readline()

	while line != '':
		split_line = line.rstrip().split('\t')
		chrom_pos = '\t'.join(split_line[:2])

		if chrom_pos in mapping:
			new_vcf.write(mapping[chrom_pos]+'\t'+'\t'.join(split_line[2:])+'\n')
		elif chrom_pos in gapmap:
			gap_vcf.write(gapmap[chrom_pos]+'\t'+'\t'.join(split_line[2:])+'\n')

		line = old_vcf.readline()

	old_vcf.close()
	new_vcf.close()
	gap_vcf.close()


# Main flow
args = parse_args()

print('Creating mapping from {0}'.format(args.delta_file))
sys.stdout.flush()

mapping, gapmap = create_mapping(args.delta_file, args.reverse)

print('Translating {0}'.format(args.VCF_file))
sys.stdout.flush()

translate_vcf(args.VCF_file, mapping, gapmap)

print('Done')

