#!/usr/bin/python

###########################################################
#
# in-core.py
#
# Determine whether a set of positions in a genome are
# within core regions
#
# Author: Connor Howington
# Date: 7/11/17
#
# Written using Python 2.6.6
#
###########################################################

import argparse
from collections import defaultdict

def parse_args():
	parser = argparse.ArgumentParser(description='Determines whether a set of positions in a VCF file are within core regions')

	parser.add_argument('core_file')

	parser.add_argument('--in_core', '-i', nargs='*', help='Lists the SNPs within the given VCF file(s) that are in core regions')
	parser.add_argument('--not_in_core', '-n', nargs='*', help='Lists the SNPs within the given VCF file(s) that are not in core regions')
	parser.add_argument('--all', '-a', nargs='*', help='Lists all SNPs in the given VCF file(s) and tells whether each is in a core region or not')

	args = parser.parse_args()

	return args


def read_cores(path):
	cores = defaultdict(lambda: [])

	with open(path, 'r') as core_file:
		line = core_file.readline().rstrip()

		while line != '':
			split_line = line.split(' ')

			chrom = split_line[0]
			start = int(split_line[1])
			end = int(split_line[2])

			cores[chrom].append([start, end, line])

			line = core_file.readline()
	return cores


def in_core(vcf_path, cores, flag):
	old_vcf_file = open(vcf_path, 'r')

	if flag == 0:
		new_vcf_path = vcf_path+'.not'
	elif flag == 1:
		new_vcf_path = vcf_path+'.core'
	else:
		new_vcf_path = vcf_path+'.all'

	new_vcf_file = open(new_vcf_path, 'w', 1)

	line = old_vcf_file.readline().rstrip()

	while line[:2] == '##':
		new_vcf_file.write(line+'\n')
		line = old_vcf_file.readline().rstrip()

	new_vcf_file.write(line+'\tIN_CORE\n')

	line = old_vcf_file.readline().rstrip()

	while line != '':
		split_line = line.split('\t')

		chrom = split_line[0]
		pos = int(split_line[1])

		found = False

		if chrom in cores:
			for core in cores[chrom]:
				if core[0] <= pos <= core[1]:
					found = True
					if flag >= 1:
						new_vcf_file.write(line+'\t'+core[2].rstrip()+'\n')
						break

		if not found and flag != 1:
			new_vcf_file.write(line+'\t.\n')

		line = old_vcf_file.readline().rstrip()

	old_vcf_file.close()
	new_vcf_file.close()



# Main flow
args = parse_args()
cores = read_cores(args.core_file)

if args.in_core is not None:
	for path in args.in_core:
		in_core(path, cores, 1)

if args.not_in_core is not None:
	for path in args.not_in_core:
		in_core(path, cores, 0)

if args.all is not None:
	for path in args.all:
		in_core(path, cores, 2)

