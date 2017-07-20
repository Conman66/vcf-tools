import sys
from collections import defaultdict

counter = defaultdict(lambda: 0)
gene_info = dict()

gene_file = open(sys.argv[1], 'r')
info_file = open(sys.argv[2], 'r')

line = gene_file.readline().rstrip()

while line != '':
	counter[line] += 1
	line = gene_file.readline().rstrip()

line = info_file.readline().rstrip()
split_line = line.split('\t')
split_line.pop(1)

print('[Our SNPs]\t'+'\t'.join(split_line))

while line != '':
	split_line = line.split('\t')
	split_line.pop(1)
	gene_info[split_line[0]] = '\t'.join(split_line)
	line = info_file.readline().rstrip()

keys = list(counter.keys())
values = list(counter.values())

while len(keys) > 0:
	value_index = values.index(max(values))
	value = values.pop(value_index)

	key = keys.pop(value_index)

	print(str(value)+'\t'+gene_info[key])

