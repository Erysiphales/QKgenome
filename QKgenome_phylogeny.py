#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
QKgenome_phylogeny.py [options] datasets coverage prefix
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
Generates a Phylip formatted alignment of polymorphic SNP sites based on output from QKgenome_conversion.py
This performs the following:
	1. import candidate gene analysis
	2. apply restriction (user specified coverage, SNPs, no Stop codon modification/indels)
	3. import FASTA files
	4. export for phylogenetic analysis with Phylip
	5. export matrix for coverage, presence of InDels, stop codon modifications, and mutations in introns for each gene x sample

Future improvements to include:
"""

## modules
import Bio
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import optparse
from optparse import OptionParser 
import sets
import string


# global variables
DNA = ['A', 'C', 'G', 'T']
IUPAC_SNP = {'R':['A', 'G'], 'Y':['C', 'T'], 'M':['A', 'C'], 'K':['G', 'T'], 'W':['A', 'T'], 'S':['C', 'G']}


## OptionParser
# import arguments and options
usage = "usage: %prog [options] datasets coverage prefix"
parser = OptionParser(usage=usage)
parser.add_option("-m", "--mask", action="store_true", dest="mask", default=False, help="Mask sequence below read coverage threshold")
parser.add_option("-s", "--synonymous", action="store_true", dest="synonymous", default=False, help="Restrict analysis to synonymous SNPs")
parser.add_option("-a", "--all", action="store_true", dest="allsites", default=False, help="Export all sites (monomorphic and polymorphic)")
(options, args) = parser.parse_args()


## exceptions
if options.mask and options.synonymous:
	print '\t' + 'Error: Commands -m or -s cannot be selected in parallel'
	exit()


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


# set user-specified coverage
coverage = float(args[1])


# read individual data sets - candidate gene analysis
print 'read individual data sets - candidate gene analysis'
gene_order = []
dataset_gene = {}
dataset_gene_with_coverage = {}
dataset_gene_coverage_indel_stop_intron = {}
genes_with_SNPs = []

for dataset in dataset_order:
	print '\t' + dataset,

	# read individual analyses from QKgenome_conversion.py
	candidate_gene_analysis_file = open(dataset + '_candidate_gene_analysis.txt', 'r')

	dataset_gene[dataset] = []
	dataset_gene_with_coverage[dataset] = []
	dataset_gene_coverage_indel_stop_intron[dataset] = {}

	truth = False

	for line in candidate_gene_analysis_file.readlines():
		line = string.replace(line, '\n', '')
		sline = string.split(line, '\t')

		if truth:
			if options.mask:
				if sline[9] == 'NA':
					dataset_gene[dataset].append(sline[0])
					genes_with_SNPs.append(sline[0])
					dataset_gene_with_coverage[dataset].append(sline[0])
			elif float(sline[5]) >= coverage:
				if sline[10] == 'No':
					if sline[9] == 'NA':
						if int(sline[6]) > 0:
							dataset_gene[dataset].append(sline[0])
							genes_with_SNPs.append(sline[0])

						dataset_gene_with_coverage[dataset].append(sline[0])

			dataset_gene_coverage_indel_stop_intron[dataset][sline[0]] = [sline[5], sline[9], sline[10], []]

			if dataset == dataset_order[0]:
				gene_order.append(sline[0])

		truth = True

	print '\t' + str(len(dataset_gene_with_coverage[dataset]))

	candidate_gene_analysis_file.close()

	intron_file = open(dataset + '_intron_junction_sequence.txt', 'r')

	for line in intron_file.readlines():
		sline = string.split(line)

		if sline[2] != sline[3]:
			dataset_gene_coverage_indel_stop_intron[dataset][sline[0]][3].append(sline[1])
	
	intron_file.close()

genes_with_coverage = dataset_gene_with_coverage[dataset_order[0]]

for index in range(1, len(dataset_order)):
	genes_with_coverage = list(sets.Set(genes_with_coverage) & sets.Set(dataset_gene_with_coverage[dataset_order[index]]))

genes_with_SNPs = list(sets.Set(genes_with_SNPs))

print
print 'Initial assessment'
print '--------------------------------------------------'
print 'Number of genes with at least 1 SNP:', len(genes_with_SNPs)
print 'Genes with sufficient coverage across datasets:', len(genes_with_coverage)
print 'Intersection of the above two sets', len(sets.Set(genes_with_SNPs) & sets.Set(genes_with_coverage))
print '--------------------------------------------------'
print
print 'dataset' + '\t' + 'shared genes with SNPs' + '\t' + 'shared genes with coverage (across) and SNPs'

for dataset in dataset_order:
	print dataset + '\t' + str(len(sets.Set(genes_with_SNPs) & sets.Set(dataset_gene[dataset]))) + '\t' + str(len(sets.Set(genes_with_coverage) & sets.Set(dataset_gene[dataset])))

print '--------------------------------------------------'
print

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


# export matrix with InDel, stop codon modifications, and intron junction mutations
coverage_file = open(args[2] + '_coverage.txt', 'w')
indel_file = open(args[2] + '_CDS_InDel.txt', 'w')
stop_codon_file = open(args[2] + '_stop_codon_modification.txt', 'w')
intron_junction_file = open(args[2] + '_intron_junctions.txt', 'w')

coverage_file.write('gene')
indel_file.write('gene')
stop_codon_file.write('gene')
intron_junction_file.write('gene')

for dataset in dataset_order:
	coverage_file.write('\t' + datasets[dataset])
	indel_file.write('\t' + datasets[dataset])
	stop_codon_file.write('\t' + datasets[dataset])
	intron_junction_file.write('\t' + datasets[dataset])

coverage_file.write('\n')
indel_file.write('\n')
stop_codon_file.write('\n')
intron_junction_file.write('\n')

for gene in gene_order:
	coverage_file.write(gene)
	indel_file.write(gene)
	stop_codon_file.write(gene)
	intron_junction_file.write(gene)

	for dataset in dataset_order:
		coverage_file.write('\t' + dataset_gene_coverage_indel_stop_intron[dataset][gene][0])
		indel_file.write('\t' + dataset_gene_coverage_indel_stop_intron[dataset][gene][1])
		stop_codon_file.write('\t' + dataset_gene_coverage_indel_stop_intron[dataset][gene][2])
		intron_junction_file.write('\t' + str(len(dataset_gene_coverage_indel_stop_intron[dataset][gene][3])))
	
	coverage_file.write('\n')
	indel_file.write('\n')
	stop_codon_file.write('\n')
	intron_junction_file.write('\n')

coverage_file.close()
indel_file.close()
stop_codon_file.close()
intron_junction_file.close()

# export sequence source information for phylogenetic analysis
print 'export sequence source information for phylogenetic analysis'
gene_position_file = open(args[2] + '_source_information.txt', 'w')

gene_position_file.write('gene' + '\t' + 'position')

if not options.mask:
	protein_position_file = open(args[2] + '_protein_variation.txt', 'w')
	protein_position_file.write('gene' + '\t' + 'position')

for dataset in dataset_order:
	gene_position_file.write('\t' + datasets[dataset])

	if not options.mask:
		protein_position_file.write('\t' + datasets[dataset])

gene_position_file.write('\n')

if not options.mask:
	protein_position_file.write('\n')

dataset_phylip_input = {}

for dataset in dataset_order:
	dataset_phylip_input[dataset] = ''

if options.allsites:
	gene_selection = list(sets.Set(genes_with_coverage))
else:
	gene_selection = list(sets.Set(genes_with_SNPs) & sets.Set(genes_with_coverage))


# current version of the script does not assess redundancy
# use code below for data set specific redundacy removal
# this will depend on the approach for naming splice model
"""
gene_redundancy = {}
gene_selection = []

for gene in list(gene_selection):
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
"""

print

gene_selection.sort()

for gene in gene_selection:
	CDS =  Seq(dataset_gene_sequence[dataset_order[0]][gene], IUPAC.ambiguous_dna)
	pep = CDS.translate(to_stop=True)

	gene_length = len(dataset_gene_sequence[dataset_order[0]][gene])
	protein_length = len(pep)

	peptides = {}

	for dataset in dataset_order:
		CDS =  Seq(dataset_gene_sequence[dataset][gene], IUPAC.ambiguous_dna)
		pep = CDS.translate(to_stop=True)

		peptides[dataset] = str(pep)

	if not options.mask:
		for aa_index in range(protein_length):
			aminoacids = []
	
			for dataset in dataset_order:
				aminoacids.append(peptides[dataset][aa_index])
	
			if len(sets.Set(aminoacids)) == 1:
				if options.synonymous:
					for base_index in range(3):
						nucleotides = []
	
						for dataset in dataset_order:
							nucleotides.append(dataset_gene_sequence[dataset][gene][aa_index * 3 + base_index])
	
						if len(sets.Set(nucleotides) - sets.Set(['N'])) > 1:
							gene_position_file.write(gene + '\t' + str(aa_index * 3 + base_index + 1))
	
							for dataset in dataset_order:
								dataset_phylip_input[dataset] += dataset_gene_sequence[dataset][gene][aa_index * 3 + base_index]
								gene_position_file.write('\t' + dataset_gene_sequence[dataset][gene][aa_index * 3 + base_index])
	
							gene_position_file.write('\n')
			else:
				protein_position_file.write(gene + '\t' + str(aa_index + 1))
	
				for dataset in dataset_order:
					protein_position_file.write('\t' + peptides[dataset][aa_index])
	
				protein_position_file.write('\n')
		
	if not options.synonymous:
		for base_index in range(gene_length):
			nucleotides = []

			for dataset in dataset_order:
				nucleotides.append(dataset_gene_sequence[dataset][gene][base_index])

			if len(sets.Set(nucleotides) - sets.Set(['N'])) > 1:
				gene_position_file.write(gene + '\t' + str(base_index + 1))

				for dataset in dataset_order:
					dataset_phylip_input[dataset] += dataset_gene_sequence[dataset][gene][base_index]
					gene_position_file.write('\t' + dataset_gene_sequence[dataset][gene][base_index])

				gene_position_file.write('\n')

if options.allsnps:
	for dataset in dataset_order:
		dataset_phylip_input[dataset] = ''

		for gene in gene_selection:
			dataset_phylip_input[dataset] += dataset_gene_sequence[dataset][gene]

gene_position_file.close()

if not options.mask:
	protein_position_file.close()

print 'number of genes in alignment:', len(gene_selection)
print 'number of SNPs in alignment:', len(dataset_phylip_input[datasets.keys()[0]])

evaluated_positions = 0

for gene in list(sets.Set(gene_selection)):
	evaluated_positions += len(dataset_gene_sequence[dataset_order[0]][gene])

print 'number of evaluated positions:', evaluated_positions, 'bp'


# export multiple sequence alignment of polymorphic sites
print 'export multiple sequence alignment of polymorphic sites'
phylip_input = open(args[2] + '_phylogeny.phy', 'w')

phylip_input.write(' ' + str(len(dataset_order)) + ' ' + str(len(dataset_phylip_input[dataset_order[0]])) + '\n')

for dataset in dataset_order:
	phylip_input.write(datasets[dataset] + ''.join([' '] * (10 - len(datasets[dataset]))) + dataset_phylip_input[dataset] + '\n')

phylip_input.close()
