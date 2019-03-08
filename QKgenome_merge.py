#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
QKgenome_merge.py [options] datasets
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
Concatenates coding sequences for several genotypes
This performs the following:
	1. import candidate gene analysis
	2. import FASTA files
	3. concatenate CDS sequences
	4. export merged CDS sequences

Future improvements to include:
"""

## modules
import optparse
from optparse import OptionParser 
import sets
import string


## OptionParser
# import arguments and options
usage = "usage: %prog [options] datasets concatenated_file [gene list]"
parser = OptionParser(usage=usage)
(options, args) = parser.parse_args()


# read in dataset identifiers
filename_file = open(args[0], 'r')

datasets = {}
dataset_order = []

for line in filename_file.readlines():
	line = string.replace(line, '\n', '')
	sline = string.split(line, '\t')

	if len(sline) > 0:
		datasets[sline[0]] = sline[1]
		dataset_order.append(sline[0])

filename_file.close()


# read in selected genes, if provided 
if len(args) > 2:
	selected_genes_file = open(args[2], 'r')
	
	selected_genes = []
	
	for line in selected_genes_file.readlines():
		sline = string.split(line)
	
		if len(sline) > 0:
			selected_genes.append(sline[0])
	
	selected_genes_file.close()


# read individual data sets - candidate gene analysis
print 'read individual data sets - candidate gene analysis'
gene_order = []

for dataset in dataset_order:
	# read individual analyses from QKgenome_conversion.py
	candidate_gene_analysis_file = open(dataset + '_candidate_gene_analysis.txt', 'r')

	truth = False

	for line in candidate_gene_analysis_file.readlines():
		line = string.replace(line, '\n', '')
		sline = string.split(line, '\t')

		if truth:
			if dataset == dataset_order[0]:
				gene_order.append(sline[0])

		truth = True

	candidate_gene_analysis_file.close()

# read individual data sets - coding sequence
print 'read individual data sets - coding sequence'
dataset_gene_sequence = {}

for dataset in dataset_order:
	cds_sequence_file = open(dataset + '_CDS.fa', 'r')

	dataset_gene_sequence[dataset] = {}

	for line in cds_sequence_file.readlines():
		if len(line) > 0:
			if line[0] == '>':
				ID = string.split(line)[0][1:]
				dataset_gene_sequence[dataset][ID] = ''
			else:
				dataset_gene_sequence[dataset][ID] += string.split(line)[0]

	cds_sequence_file.close()


# current version of the script does not assess redundancy
# use code below for data set specific redundacy removal
# this will depend on the approach for naming splice model
gene_redundancy = {}
gene_selection = []

for gene in list(sets.Set(gene_order)):
	if string.split(gene, '.')[0] not in gene_redundancy.keys():
		gene_redundancy[string.split(gene, '.')[0]] = []
	
	gene_redundancy[string.split(gene, '.')[0]].append(gene)

for gene in gene_redundancy.keys():
	if len(gene_redundancy[gene]) == 1:
		gene_selection.append(gene_redundancy[gene][0])
	elif len(gene_redundancy[gene]) > 1:
		selected = False

		for splice_model in gene_redundancy[gene]:
			if '.1.v3.1' in splice_model:
				gene_selection.append(splice_model)
				selected = True

		if not selected:
			for splice_model in gene_redundancy[gene]:
				if '.2.v3.1' in splice_model:
					gene_selection.append(splice_model)
					selected = True

		if not selected:
			print '\tFreak out!'

		print gene, len(gene_redundancy[gene]), gene_redundancy[gene] 

gene_selection.sort()

# concatenate and export CDS sequence
concatenated_file = open(args[1], 'w')

for dataset in dataset_order:
	dataset_concat_sequence = ''

	if len(args) > 2:
		for gene in selected_genes:
			dataset_concat_sequence += dataset_gene_sequence[dataset][gene]
	else:
		for gene in gene_selection:
			dataset_concat_sequence += dataset_gene_sequence[dataset][gene]
	
	concatenated_file.write('>' + datasets[dataset] + '\n')
	concatenated_file.write(dataset_concat_sequence + '\n')

concatenated_file.close()
