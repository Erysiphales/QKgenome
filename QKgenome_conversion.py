#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog [options] coverage_threshold frequency_threshold FASTA GFF3 pileup2snp pileup2indel genomecov annotation prefix [expression files]
Expression files from featureCounts are optional.
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
Convert a genome sequence and annotations based on read alignment
This performs the following:
	1. <insert protocol>

Improvements from existing script set (9 February 2017):
	1. Fixed InDel sign with respect to reference
	2. Initialize with prefix, skip individual file output
	3. Remove "Interupted ORF column", more descriptive terms
	4. Add analysis of intron splice junctions, specifically looking for SNPs at these positions
	5. Export final file based on start position of ORFs (skips ordering step in Excel)
	6. Add ability to plot figures that show
		a. Variant frequency (with coverage threshold)
		b. Average coverage across entire data set (gene space only)
		c. Frequency of different splice sites in region
		d. Candidate gene analysis file, well structured (use NA), export relevant images

Future improvements to include:
	1. Incorporate dN/dS analysis?
	2. Add ability to investigate heterokaryotic SNPs (such as in dikaryotic rusts)
	3. Export statistics and information about alignment, with call information as well (for book keeping purposes, with date/time)
		a. Document time - command is time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
		b. Number of genes evaluated
		c. ???
"""

## modules
import Bio
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import commands
import time
import optparse
from optparse import OptionParser 
import sets
import string

"""
# incorporate later to have dN/dS analyses

import synonymous_calc_MM1
from synonymous_calc_MM1 import *
"""

## global variables
DNA = ['A', 'C', 'G', 'T']
cmp_DNA = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'M':'K', 'R':'Y', 'W':'S', 'S':'W', 'Y':'R', 'K':'M', 'V':'B', 'H':'D', 'D':'H', 'B':'V', 'X':'X', 'N':'N'}
IUPAC_convert = {'AG':'R', 'GA':'R', 'CT':'Y', 'TC':'Y', 'CA':'M', 'AC':'M', 'TG':'K', 'GT':'K', 'TA':'W', 'AT':'W', 'CG':'S', 'GC':'S'}
IUPAC_SNP = {'R':['A', 'G'], 'Y':['C', 'T'], 'M':['A', 'C'], 'K':['G', 'T'], 'W':['A', 'T'], 'S':['C', 'G']}


## functions
# reverse complement
# input : DNA sequence
# output : reverse complement of said DNA sequence
def reverse_complement(orig_DNA):
	rev_comp_DNA = ''
	for index in range(len(orig_DNA)):
		rev_comp_DNA += cmp_DNA[orig_DNA[len(orig_DNA)-index-1]]
	return rev_comp_DNA

# import_SNPs
# reads through a VarScan output file and only retains SNPs that are above the read count and percent threshold
# contig -> position of SNP -> [reference allele, alternate allele, percent variant frequency]
# modified to work with differing output from VarScan (pileup versus mpileup)
def import_SNPs(varscan, read_threshold, percent_threshold, prefix):
	contig_position_allele = {}

	varscan_file = open(varscan, 'r')
	SNP_frequency_file = open(prefix + '_SNP_frequency.txt', 'w')
	SNP_frequency_file.write('coverage' + '\t' + 'frequency' + '\n')
	
	truth = False

	for line in varscan_file.readlines():
		sline = string.split(line)

		if truth:
			if '%' in sline[4]:
				if int(string.split(sline[4], ':')[1]) >= read_threshold:
					SNP_frequency_file.write(string.split(sline[4], ':')[1] + '\t' + string.replace(string.split(sline[4], ':')[4], '%', '') + '\n')

					if float(string.replace(string.split(sline[4], ':')[4], '%', '')) >= percent_threshold:
						if sline[0] not in contig_position_allele.keys():
							contig_position_allele[sline[0]] = {}

						contig_position_allele[sline[0]][int(sline[1])] = [sline[2], sline[3], float(string.replace(string.split(sline[4], ':')[4], '%', ''))]
			else:
				if (int(sline[4]) + int(sline[5])) >= read_threshold:
					SNP_frequency_file.write(str(int(sline[4]) + int(sline[5])) + '\t' + string.replace(sline[6], '%', '') + '\n')

					if float(string.replace(sline[6], '%', '')) >= percent_threshold:
						if sline[0] not in contig_position_allele.keys():
							contig_position_allele[sline[0]] = {}

						contig_position_allele[sline[0]][int(sline[1])] = [sline[2], sline[18], float(string.replace(sline[6], '%', ''))]

		truth = True

	varscan_file.close()
	SNP_frequency_file.close()

	return contig_position_allele

# import_heterozygous_SNPs
# reads through a VarScan output file and only retains SNPs that are above the read count and between the lower and upper percent thresholds
# contig -> position of SNP -> [reference allele, alternate allele, percent variant frequency]
# modified to work with differing output from VarScan (pileup versus mpileup)
def import_heterozygous_SNPs(varscan, read_threshold, lower_threshold, higher_threshold, percent_threshold, prefix):
	contig_position_allele = {}

	varscan_file = open(varscan, 'r')
	SNP_frequency_file = open(prefix + '_heterozygous_SNP_frequency.txt', 'w')
	SNP_frequency_file.write('coverage' + '\t' + 'frequency' + '\n')
	
	truth = False

	for line in varscan_file.readlines():
		sline = string.split(line)

		if truth:
			if '%' in sline[4]:
				if int(string.split(sline[4], ':')[1]) >= read_threshold:
					SNP_frequency_file.write(string.split(sline[4], ':')[1] + '\t' + string.replace(string.split(sline[4], ':')[4], '%', '') + '\n')

					if float(string.replace(string.split(sline[4], ':')[4], '%', '')) >= lower_threshold:
						if float(string.replace(string.split(sline[4], ':')[4], '%', '')) < higher_threshold:
							if sline[0] not in contig_position_allele.keys():
								contig_position_allele[sline[0]] = {}

							if int(sline[1]) not in contig_position_allele[sline[0]].keys():
								if sline[3] in IUPAC_SNP.keys():
									variant_allele = list(sets.Set(IUPAC_SNP[sline[3]]) - sets.Set(sline[2]))[0]
									contig_position_allele[sline[0]][int(sline[1])] = [sline[2], variant_allele, float(string.replace(string.split(sline[4], ':')[4], '%', ''))]
								else:
									contig_position_allele[sline[0]][int(sline[1])] = [sline[2], sline[3], float(string.replace(string.split(sline[4], ':')[4], '%', ''))]
							else:
								variant_allele = list(sets.Set(IUPAC_SNP[sline[3]]) - sets.Set(sline[2]))[0]
								true_allele = IUPAC_convert[contig_position_allele[sline[0]][int(sline[1])][1] + variant_allele]
								contig_position_allele[sline[0]][int(sline[1])][1] = true_allele
			else:
				if (int(sline[4]) + int(sline[5])) >= read_threshold:
					SNP_frequency_file.write(str(int(sline[4]) + int(sline[5])) + '\t' + string.replace(sline[6], '%', '') + '\n')

					if float(string.replace(sline[6], '%', '')) >= lower_threshold:
						if float(string.replace(sline[6], '%', '')) < higher_threshold:
							if sline[0] not in contig_position_allele.keys():
								contig_position_allele[sline[0]] = {}

							if int(sline[1]) not in contig_position_allele[sline[0]].keys():
								if sline[3] in IUPAC_SNP.keys():
									variant_allele = list(sets.Set(IUPAC_SNP[sline[3]]) - sets.Set(sline[2]))[0]
									contig_position_allele[sline[0]][int(sline[1])] = [sline[2], variant_allele, float(string.replace(sline[6], '%', ''))]
								else:
									contig_position_allele[sline[0]][int(sline[1])] = [sline[2], sline[3], float(string.replace(sline[6], '%', ''))]
							else:
								variant_allele = list(sets.Set(IUPAC_SNP[sline[3]]) - sets.Set(sline[2]))[0]
								true_allele = IUPAC_convert[contig_position_allele[sline[0]][int(sline[1])][1] + variant_allele]
								contig_position_allele[sline[0]][int(sline[1])][1] = true_allele

		truth = True

	varscan_file.close()
	SNP_frequency_file.close()

	return contig_position_allele

# import_indels
# reads through a VarScan output file and only retains indels that are above the read count and percent threshold
# contig -> position of indel -> [reference allele, insertion/deletion, alternate allele, percent variant frequency]
# modified to work with differing output from VarScan (pileup versus mpileup)
def import_indels(varscan, read_threshold, percent_threshold, prefix):
	contig_position_allele = {}

	varscan_file = open(varscan, 'r')
	indel_frequency_file = open(prefix + '_indel_frequency.txt', 'w')
	indel_frequency_file.write('coverage' + '\t' + 'frequency' + '\n')

	truth = False

	for line in varscan_file.readlines():
		sline = string.split(line)

		if truth:
			if '%' in sline[4]:
				if int(string.split(sline[4], ':')[1]) >= read_threshold:
					indel_frequency_file.write(string.split(sline[4], ':')[1] + '\t' + string.replace(string.split(sline[4], ':')[4], '%', '') + '\n')

					if float(string.replace(string.split(sline[4], ':')[4], '%', '')) >= percent_threshold:
						if sline[0] not in contig_position_allele.keys():
							contig_position_allele[sline[0]] = {}

						contig_position_allele[sline[0]][int(sline[1])] = [sline[2], sline[3][0], sline[3][1:], float(string.replace(string.split(sline[4], ':')[4], '%', ''))]
			else:
				if (int(sline[4]) + int(sline[5])) >= read_threshold:
					indel_frequency_file.write(str(int(sline[4]) + int(sline[5])) + '\t' + string.replace(sline[6], '%', '') + '\n')

					if float(string.replace(sline[6], '%', '')) >= percent_threshold:
						if sline[0] not in contig_position_allele.keys():
							contig_position_allele[sline[0]] = {}

						contig_position_allele[sline[0]][int(sline[1])] = [sline[2], sline[18][0], sline[18][1:], float(string.replace(sline[6], '%', ''))]

		truth = True

	varscan_file.close()
	indel_frequency_file.close()

	return contig_position_allele

# str_distance
# determine the distance between two strings (useful for nucleotide and protein)
# if lengths of strings are different or if either string is of length 0, return -1
# string1, string2 -> integer
def str_distance(string1, string2):
	if len(string1) != len(string2):
		distance = -1
	elif len(string1) == 0:
		distance = -1
	elif len(string2) == 0:
		distance = -1
	else:
		distance = 0

		for index in range(len(string1)):
			if string1[index] != string2[index]:
				distance += 1
	
	return distance


## OptionParser
# import arguments and options
usage = "usage: %prog [options] coverage_threshold frequency_threshold FASTA GFF3 pileup2snp pileup2indel genomecov prefix [expression files]"
parser = OptionParser(usage=usage)
parser.add_option("-a", "--annotation", action="store", type="string", dest="annotation", default='', help="Annotation file (tab-limited)")
parser.add_option("-d", "--heterozygous", action="store_true", dest="het", default=False, help="Evaluate heterozygous/hemizygous/dikaryotic SNPs and InDels")
parser.add_option("-m", "--mask", action="store_true", dest="mask", default=False, help="Mask sequence below read coverage threshold")
parser.add_option("-l", "--lower", action="store", type="float", dest="lower", default=-1, help="Lower threshold for heterozygous SNPs")
parser.add_option("-u", "--upper", action="store", type="float", dest="upper", default=-1, help="Upper threshold for heterozygous SNPs")
(options, args) = parser.parse_args()


# coverage and frequency parameters
coverage_threshold = int(args[0])
variant_frequency_threshold = float(args[1])
prefix = args[7]
heterozygous_analysis = False

summary_file = open(prefix + '_analysis_' + time.strftime("%Y%m%d_%Hh%Mm%Ss", time.gmtime()) + '.txt', 'w')
summary_file.write('python %prog')

if options.mask:
	summary_file.write(' -m')
if len(options.annotation) > 0:
	summary_file.write(' -a ' + options.annotation)
if options.lower > 0:
	if options.lower < 100:
		if options.upper > 0:
			if options.upper < 100:
				summary_file.write(' -l ' + str(options.lower) + ' -u ' + str(options.upper))
				heterozygous_analysis = True


summary_file.write(' ' + ' '.join(args) + '\n')
summary_file.write('Run start time: ' + time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + '\n')

# read in reference genome, determine size of all contigs
print 'FASTA'
summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'fasta' + '\n')
fasta = open(args[2], 'r')

ID_sequence = {}
ID_sequence_converted = {}
contig_adjustment_base_pair = {}

for line in fasta.readlines():
	if len(line) > 0:
		if line[0] == '>':
			ID = string.split(line)[0][1:]
			ID_sequence[ID] = ''
			contig_adjustment_base_pair[ID] = []
		else:
			ID_sequence[ID] += string.split(line)[0]

fasta.close()

for ID in ID_sequence.keys():
	ID_sequence_converted[ID] = ID_sequence[ID]

# import SNPs
print 'SNPs'
summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'SNPs' + '\n')
contig_position_allele_SNPs = import_SNPs(args[4], coverage_threshold, variant_frequency_threshold, prefix)

if heterozygous_analysis:
	contig_position_allele_heterozygous_SNPs = import_heterozygous_SNPs(args[4], coverage_threshold, options.lower, options.upper, variant_frequency_threshold, prefix)

# import indels
print 'InDels'
summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'InDels' + '\n')
contig_position_allele_indels = import_indels(args[5], coverage_threshold, variant_frequency_threshold, prefix)

# processing genomic sequence
print 'Convert SNPs'
summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'Convert SNPs' + '\n')
for contig in contig_position_allele_SNPs.keys():
	for position in contig_position_allele_SNPs[contig].keys():
		ID_sequence_converted[contig] = ID_sequence_converted[contig][:(position - 1)] + contig_position_allele_SNPs[contig][position][1] + ID_sequence_converted[contig][position:]

if heterozygous_analysis:
	for contig in contig_position_allele_heterozygous_SNPs.keys():
		for position in contig_position_allele_heterozygous_SNPs[contig].keys():
			ID_sequence_converted[contig] = ID_sequence_converted[contig][:(position - 1)] + contig_position_allele_heterozygous_SNPs[contig][position][1] + ID_sequence_converted[contig][position:]

print 'Convert InDels'
summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'Convert InDels' + '\n')
for contig in contig_position_allele_indels.keys():
	base_pair_offset = 0

	for position in range(len(ID_sequence_converted[contig])):
		contig_adjustment_base_pair[contig].append(base_pair_offset)

		if position in contig_position_allele_indels[contig].keys():
			if contig_position_allele_indels[contig][position][1] == '+':
				ID_sequence_converted[contig] = ID_sequence_converted[contig][:(position + base_pair_offset)] + contig_position_allele_indels[contig][position][2] + ID_sequence_converted[contig][(position + base_pair_offset):]

				base_pair_offset += len(contig_position_allele_indels[contig][position][2])
			elif contig_position_allele_indels[contig][position][1] == '-':
				ID_sequence_converted[contig] = ID_sequence_converted[contig][:(position + base_pair_offset)] + ID_sequence_converted[contig][(position + base_pair_offset + len(contig_position_allele_indels[contig][position][2])):]

				base_pair_offset -= len(contig_position_allele_indels[contig][position][2])
			else:
				print 'Here be dragons...'

base_pair_offset = 0

for contig in list(sets.Set(ID_sequence_converted.keys()) - sets.Set(contig_position_allele_indels.keys())):
	for position in range(len(ID_sequence_converted[contig])):
		contig_adjustment_base_pair[contig].append(base_pair_offset)

print 'Export Modified Sequence'
summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'Export Modified Sequence' + '\n')

converted_sequence_file = open(prefix + '.fa', 'w')

for ID in ID_sequence_converted.keys():
	converted_sequence_file.write('>' + ID + '\n')
	converted_sequence_file.write(ID_sequence_converted[ID] + '\n')

converted_sequence_file.close()


# read in GFF3, write adjusted GFF3
print 'GFF3'
summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'GFF3' + '\n')

contig_geneID_transcript_start_stop_strand = {}
contig_geneID_transcript_start_stop_strand_converted = {}
geneID_exons = {}
geneID_exons_converted = {}
geneID_CDS_exons = {}
geneID_CDS_exons_converted = {}
contig_gene_clusters = {}
geneID_CDS_start_stop = {}
geneID_CDS_start_stop_converted = {}

# initialize contigs
GTF_file = open(args[3], 'r')

for line in GTF_file.readlines():
	if len(line) > 0:
		if line[0] != '#':
			sline = string.replace(line, '\n', '')
			sline = string.split(line, '\t')

			contig_geneID_transcript_start_stop_strand[sline[0]] = {}
			contig_geneID_transcript_start_stop_strand_converted[sline[0]] = {}
	
GTF_file.close()

# reads information
GTF_file = open(args[3], 'r')
converted_GTF_file = open(prefix + '.gff3', 'w')
geneID = ''

for line in GTF_file.readlines():
	if len(line) > 0:
		if line[0] != '#':
			line = string.replace(line, '\n', '')
			sline = string.split(line, '\t')

			contig = sline[0]
		
			if sline[2] == 'mRNA':
				geneID = string.split(sline[8], ';')[0][3:]

				contig_geneID_transcript_start_stop_strand[sline[0]][geneID] = [int(sline[3]) - 1, int(sline[4]) - 1, sline[6]]
				contig_geneID_transcript_start_stop_strand_converted[sline[0]][geneID] = [int(sline[3]) - 1 + contig_adjustment_base_pair[contig][int(sline[3]) - 1], int(sline[4]) - 1 + contig_adjustment_base_pair[contig][int(sline[4]) - 1], sline[6]]

				geneID_exons[geneID] = []
				geneID_exons_converted[geneID] = []

				geneID_CDS_exons[geneID] = []
				geneID_CDS_exons_converted[geneID] = []

			if sline[2] == 'CDS':
				geneID = string.split(sline[8], ';')[1][7:]

				geneID_exons[geneID].append([int(sline[3]) - 1, int(sline[4]) - 1])
				geneID_exons_converted[geneID].append([int(sline[3]) - 1 + contig_adjustment_base_pair[contig][int(sline[3]) - 1], int(sline[4]) - 1 + contig_adjustment_base_pair[contig][int(sline[4]) - 1]])

			if sline[2] == 'CDS':
				geneID = string.split(sline[8], ';')[1][7:]

				geneID_CDS_exons[geneID].append([int(sline[3]) - 1, int(sline[4]) - 1])
				geneID_CDS_exons_converted[geneID].append([int(sline[3]) - 1 + contig_adjustment_base_pair[contig][int(sline[3]) - 1], int(sline[4]) - 1 + contig_adjustment_base_pair[contig][int(sline[4]) - 1]])
		
			if len(geneID) > 0:
				geneID_CDS_start_stop[geneID] = {}
				geneID_CDS_start_stop_converted[geneID] = {}
			
			converted_GTF_file.write(sline[0] + '\t' + sline[1] + '\t' + sline[2] + '\t' + str(int(sline[3]) + contig_adjustment_base_pair[contig][int(sline[3]) - 1]) + '\t' + str(int(sline[4]) + contig_adjustment_base_pair[contig][int(sline[4]) - 1]) + '\t' + sline[5] + '\t' + sline[6] + '\t' + sline[7] + '\t' + sline[8] + '\n')
		else:
			converted_GTF_file.write(line)

GTF_file.close()
converted_GTF_file.close()

# import annotations
if len(options.annotation) > 0:
	print 'Annotations'
	summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'Annotations' + '\n')

	annotation_file = open(options.annotation, 'r')
	
	gene_annotation = {}
	
	for line in annotation_file.readlines():
		line = string.replace(line, '\n', '')
		sline = string.split(line, '\t')
	
		gene_annotation[sline[0]] = sline[1]
	
	annotation_file.close()

# import genomecov to estimate coverage over length of coding sequences
print 'Genomecov'
summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'Genomecov' + '\n')

genomecov_file = open(args[6], 'r')

# initialize dictionary contig_evaluated_positions with positions of genes
contig_evaluated_positions = {}

for contig in contig_geneID_transcript_start_stop_strand.keys():
	contig_evaluated_positions[contig] = []

	for geneID in contig_geneID_transcript_start_stop_strand[contig].keys():
		exon_order = range(len(geneID_CDS_exons[geneID]))

		for exon_index in exon_order:
			for position_index in range(geneID_CDS_exons[geneID][exon_index][0] + 1, geneID_CDS_exons[geneID][exon_index][1] + 2):
				contig_evaluated_positions[contig].append(position_index)

print 'Evaluated Positions'
summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'Evaluated Positions' + '\n')

contig_position_index = {}
contig_position_coverage = {}

for contig in contig_evaluated_positions.keys():
	contig_evaluated_positions[contig] = list(sets.Set(contig_evaluated_positions[contig]))
	contig_evaluated_positions[contig].sort()
	contig_position_index[contig] = 0
	contig_position_coverage[contig] = {}

	for position_index in range(len(ID_sequence[contig])):
		contig_position_coverage[contig][position_index] = 0

for line in genomecov_file.readlines():
	sline = string.split(line)

	if sline[0] in contig_evaluated_positions.keys():
		if contig_position_index[sline[0]] < len(contig_evaluated_positions[sline[0]]):
			while int(sline[1]) > contig_evaluated_positions[sline[0]][contig_position_index[sline[0]]]:
				contig_position_index[sline[0]] += 1

				if contig_position_index[sline[0]] > len(contig_evaluated_positions[sline[0]]):
					break

		if contig_position_index[sline[0]] < len(contig_evaluated_positions[sline[0]]):
			if int(sline[1]) == contig_evaluated_positions[sline[0]][contig_position_index[sline[0]]]:
				bases = int(float(sline[2])) 

				# save coverage, position in Python start (0) not GFF (1), positions with no coverage will not be present as a key
				contig_position_coverage[sline[0]][int(sline[1]) - 1] = bases

				contig_position_index[sline[0]] += 1

genomecov_file.close()


# generate a masked sequence based on coverage
if options.mask:
	print 'Mask Sequence'
	summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'Mask Sequence' + '\n')

	ID_sequence_masked = {}

	# initialize sequence
	for ID in ID_sequence.keys():
		ID_sequence_masked[ID] = ID_sequence[ID]

	# mask sequence based on coverage
	for contig in contig_position_coverage.keys():
		for position in contig_position_coverage[contig].keys():
			if contig_position_coverage[contig][position] < coverage_threshold:
				ID_sequence_masked[contig] = ID_sequence_masked[contig][:(position - 1)] + 'N' + ID_sequence_masked[contig][position:]

	# incorporate indels
	for contig in contig_position_allele_indels.keys():
		base_pair_offset = 0

		for position in range(len(ID_sequence_masked[contig])):
			if position in contig_position_allele_indels[contig].keys():
				if contig_position_allele_indels[contig][position][1] == '+':
					ID_sequence_masked[contig] = ID_sequence_masked[contig][:(position + base_pair_offset)] + contig_position_allele_indels[contig][position][2] + ID_sequence_masked[contig][(position + base_pair_offset):]

					base_pair_offset += len(contig_position_allele_indels[contig][position][2])
				elif contig_position_allele_indels[contig][position][1] == '-':
					ID_sequence_masked[contig] = ID_sequence_masked[contig][:(position + base_pair_offset)] + ID_sequence_masked[contig][(position + base_pair_offset + len(contig_position_allele_indels[contig][position][2])):]

					base_pair_offset -= len(contig_position_allele_indels[contig][position][2])
				else:
					print 'Here be dragons...'

dataset_gene_expression = {}

if len(args) > 9:
	print 'Expression'
	summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'Expression' + '\n')

	for expression_filename in args[9:]:
		expression_file = open(expression_filename, 'r')

		dataset_gene_expression[string.split(expression_filename, '_')[1]] = {}

		for line in expression_file.readlines():
			sline = string.split(line)

			geneID = string.replace(sline[0], '.v3.1', '')

			dataset_gene_expression[string.split(expression_filename, '_')[1]][geneID] = sline[6]
			
		expression_file.close()

experimental_datasets = dataset_gene_expression.keys()
experimental_datasets.sort()

# export transcript files for transdecoder analysis
print 'Export'
summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'export' + '\n')

#ORF_file = open(prefix + '_ORF.fa', 'w')
ORF_converted_file = open(prefix + '_CDS.fa', 'w')
gene_model_analysis_file = open(prefix + '_candidate_gene_analysis.txt', 'w')
intron_junction_file = open(prefix + '_intron_junction_sequence.txt', 'w')

#gene_model_analysis_file.write('Gene' + '\t' + 'ORF Start' + '\t' + 'ORF Stop' + '\t' + 'Strand' + '\t' + 'mRNA length' + '\t' + 'Alignment coverage' + '\t' + 'SNPs' + '\t' + 'Protein length' + '\t' + 'Amino acid differences' + '\t' + 'InDels' + '\t' + 'Interrupted ORF' + '\t' + 'Annotation' + '\t' + 'Pfam domain')
if len(options.annotation) > 0:
	gene_model_analysis_file.write('Gene' + '\t' + 'ORFStart' + '\t' + 'ORFStop' + '\t' + 'Strand' + '\t' + 'mRNAlength' + '\t' + 'Coverage' + '\t' + 'SNPs' + '\t' + 'Proteinlength' + '\t' + 'Aminoaciddifferences' + '\t' + 'InDels' + '\t' + 'InterruptedORF' + '\t' + 'Annotation')
else:
	gene_model_analysis_file.write('Gene' + '\t' + 'ORFStart' + '\t' + 'ORFStop' + '\t' + 'Strand' + '\t' + 'mRNAlength' + '\t' + 'Coverage' + '\t' + 'SNPs' + '\t' + 'Proteinlength' + '\t' + 'Aminoaciddifferences' + '\t' + 'InDels' + '\t' + 'InterruptedORF')

for index in range(len(experimental_datasets)):
	gene_model_analysis_file.write('\t' + experimental_datasets[index])

gene_model_analysis_file.write('\n')

total_coverage = [0, 0]

for contig in contig_geneID_transcript_start_stop_strand.keys():
	position_gene = {}

	for geneID in contig_geneID_transcript_start_stop_strand[contig].keys():
		if geneID_CDS_exons[geneID][0][0] not in position_gene.keys():
			position_gene[geneID_CDS_exons[geneID][0][0]] = []

		position_gene[geneID_CDS_exons[geneID][0][0]].append(geneID)
	
	positions = position_gene.keys()
	positions.sort()

	for position in positions:
		geneIDs = position_gene[position]
		geneIDs.sort()

		for geneID in geneIDs:
			gene_model_analysis_file.write(geneID)
	
			exon_sequence = ''
			exon_sequence_converted = ''
	
			coding_sequence = ''
			coding_sequence_converted = ''
	
			# export index of start and stop of ORF
			gene_model_analysis_file.write('\t' + str(geneID_CDS_exons[geneID][0][0] + 1))
			gene_model_analysis_file.write('\t' + str(geneID_CDS_exons[geneID][len(geneID_CDS_exons[geneID]) - 1][1] + 1))
	
			# export strand
			gene_model_analysis_file.write('\t' + contig_geneID_transcript_start_stop_strand[contig][geneID][2])
	
	
			# export mRNA length
			if contig_geneID_transcript_start_stop_strand[contig][geneID][2] == '+':
				exon_order = range(len(geneID_exons[geneID]))
			elif contig_geneID_transcript_start_stop_strand[contig][geneID][2] == '-':
				exon_order = range(len(geneID_exons[geneID]))
				exon_order.reverse()
	
			for exon_index in exon_order:
				exon_sequence += ID_sequence[contig][geneID_exons[geneID][exon_index][0]:geneID_exons[geneID][exon_index][1] + 1]
				exon_sequence_converted += ID_sequence_converted[contig][geneID_exons_converted[geneID][exon_index][0]:geneID_exons_converted[geneID][exon_index][1] + 1]
	
				# export intron junction
				if contig_geneID_transcript_start_stop_strand[contig][geneID][2] == '+':
					if exon_index < (len(exon_order) - 1):
						intron_junction_file.write(geneID + '\t' + str(exon_index) + '-' + str(exon_index + 1) + '\t' + ID_sequence[contig][(geneID_exons[geneID][exon_index][1] + 1):(geneID_exons[geneID][exon_index][1] + 3)] + ID_sequence[contig][(geneID_exons[geneID][exon_index + 1][0] - 2):geneID_exons[geneID][exon_index + 1][0]])
	
						if options.mask:
							intron_junction_file.write('\t' + ID_sequence_masked[contig][(geneID_exons_converted[geneID][exon_index][1] + 1):(geneID_exons_converted[geneID][exon_index][1] + 3)] + ID_sequence_masked[contig][(geneID_exons_converted[geneID][exon_index + 1][0] - 2):geneID_exons_converted[geneID][exon_index + 1][0]] + '\n')
						else:
							intron_junction_file.write('\t' + ID_sequence_converted[contig][(geneID_exons_converted[geneID][exon_index][1] + 1):(geneID_exons_converted[geneID][exon_index][1] + 3)] + ID_sequence_converted[contig][(geneID_exons_converted[geneID][exon_index + 1][0] - 2):geneID_exons_converted[geneID][exon_index + 1][0]] + '\n')
				if contig_geneID_transcript_start_stop_strand[contig][geneID][2] == '-':
					if exon_index > 0:
						intron_junction_file.write(geneID + '\t' + str(len(geneID_exons[geneID]) - exon_index) + '-' + str(len(geneID_exons[geneID]) - exon_index + 1) + '\t' + reverse_complement(ID_sequence[contig][(geneID_exons[geneID][exon_index - 1][1] + 1):(geneID_exons[geneID][exon_index - 1][1] + 3)] + ID_sequence[contig][(geneID_exons[geneID][exon_index][0] - 2):geneID_exons[geneID][exon_index][0]]))	
						if options.mask:
							intron_junction_file.write('\t' + reverse_complement(ID_sequence_masked[contig][(geneID_exons_converted[geneID][exon_index - 1][1] + 1):(geneID_exons_converted[geneID][exon_index - 1][1] + 3)] + ID_sequence_masked[contig][(geneID_exons_converted[geneID][exon_index][0] - 2):geneID_exons_converted[geneID][exon_index][0]]) + '\n')	
						else:
							intron_junction_file.write('\t' + reverse_complement(ID_sequence_converted[contig][(geneID_exons_converted[geneID][exon_index - 1][1] + 1):(geneID_exons_converted[geneID][exon_index - 1][1] + 3)] + ID_sequence_converted[contig][(geneID_exons_converted[geneID][exon_index][0] - 2):geneID_exons_converted[geneID][exon_index][0]]) + '\n')	
	
			gene_model_analysis_file.write('\t' + str(len(exon_sequence_converted)))
	
			# determine coverage of reads (based on threshold)
			coverage = []
	
			# correct code below to just change exon order, otherwise all code is identical and not needed
			if contig_geneID_transcript_start_stop_strand[contig][geneID][2] == '+':
				exon_order = range(len(geneID_CDS_exons[geneID]))
	
				for exon_index in exon_order:
					coding_sequence += ID_sequence[contig][geneID_CDS_exons[geneID][exon_index][0]:geneID_CDS_exons[geneID][exon_index][1] + 1]
	
					for position_index in range(geneID_CDS_exons[geneID][exon_index][0], geneID_CDS_exons[geneID][exon_index][1] + 1):
						if contig_position_coverage[contig][position_index] >= coverage_threshold:
							coverage.append(1)
							total_coverage[1] += 1
						else:
							coverage.append(0)
							total_coverage[0] += 1
				
				for exon_index in exon_order:
					if options.mask:
						coding_sequence_converted += ID_sequence_masked[contig][geneID_CDS_exons_converted[geneID][exon_index][0]:geneID_CDS_exons_converted[geneID][exon_index][1] + 1]
					else:
						coding_sequence_converted += ID_sequence_converted[contig][geneID_CDS_exons_converted[geneID][exon_index][0]:geneID_CDS_exons_converted[geneID][exon_index][1] + 1]
	
				#ORF_file.write('>' + geneID + '_source\n')
				#ORF_file.write(coding_sequence + '\n')
				ORF_converted_file.write('>' + geneID + '\n')
				ORF_converted_file.write(coding_sequence_converted + '\n')
	
				CDS =  Seq(coding_sequence, IUPAC.unambiguous_dna)
				CDSc =  Seq(coding_sequence_converted, IUPAC.unambiguous_dna)
	
				if len(sets.Set(['N', 'n']) & sets.Set(coding_sequence)) > 0:
					pep = ''
				elif (len(coding_sequence) % 3) == 0:
					pep = CDS.translate(to_stop=True)
				else:
					pep = ''
	
				if len(sets.Set(['N', 'n']) & sets.Set(coding_sequence_converted)) > 0:
					pep = ''
				elif (len(coding_sequence_converted) % 3) == 0:
					pepc = CDSc.translate(to_stop=True)
				else:
					pepc = ''
	
				#if str(pep) != str(pepc):
				#	print geneID + '\t' + 'different' + '\t' + str(len(coding_sequence) % 3) + '\t' + str(len(coding_sequence_converted) % 3)
				#else:
				#	print geneID + '\t' +  'same' + '\t' + str(len(coding_sequence) % 3) + '\t' + str(len(coding_sequence_converted) % 3)
	
			if contig_geneID_transcript_start_stop_strand[contig][geneID][2] == '-':
				exon_order = range(len(geneID_CDS_exons[geneID]))
				exon_order.reverse()
	
				for exon_index in exon_order:
					coding_sequence += ID_sequence[contig][geneID_CDS_exons[geneID][exon_index][0]:geneID_CDS_exons[geneID][exon_index][1] + 1]
	
					for position_index in range(geneID_CDS_exons[geneID][exon_index][0], geneID_CDS_exons[geneID][exon_index][1] + 1):
						if contig_position_coverage[contig][position_index] >= coverage_threshold:
							coverage.append(1)
							total_coverage[1] += 1
						else:
							coverage.append(0)
							total_coverage[0] += 1
				
				for exon_index in exon_order:
					if options.mask:
						coding_sequence_converted += ID_sequence_masked[contig][geneID_CDS_exons_converted[geneID][exon_index][0]:geneID_CDS_exons_converted[geneID][exon_index][1] + 1]
					else:
						coding_sequence_converted += ID_sequence_converted[contig][geneID_CDS_exons_converted[geneID][exon_index][0]:geneID_CDS_exons_converted[geneID][exon_index][1] + 1]
	
				#ORF_file.write('>' + geneID + '\n')
				#ORF_file.write(reverse_complement(coding_sequence) + '\n')
				ORF_converted_file.write('>' + geneID + '\n')
				ORF_converted_file.write(reverse_complement(coding_sequence_converted) + '\n')
	
				CDS =  Seq(reverse_complement(coding_sequence), IUPAC.unambiguous_dna)
				CDSc =  Seq(reverse_complement(coding_sequence_converted), IUPAC.unambiguous_dna)
	
				if (len(coding_sequence) % 3) == 0:
					pep = CDS.translate(to_stop=True)
				else:
					pep = ''
	
				pepc = ''
	
				if options.mask:
					pepc = ''
				elif (len(coding_sequence_converted) % 3) == 0:
					pepc = CDSc.translate(to_stop=True)
				else:
					pepc = ''
	
				#if str(pep) != str(pepc):
				#	print geneID + '\t' + 'different' + '\t' + str(len(coding_sequence) % 3) + '\t' + str(len(coding_sequence_converted) % 3)
				#else:
				#	print geneID + '\t' + 'same' + '\t' + str(len(coding_sequence) % 3) + '\t' + str(len(coding_sequence_converted) % 3)
	
			# export coverage
			gene_model_analysis_file.write('\t' + str(float(sum(coverage)) / len(coverage)))
	
			# if an InDel does not exist, determine the number of SNPs
			DNA_distance = str_distance(coding_sequence, coding_sequence_converted)
	
			if (len(coding_sequence) % 3) > 0:
				print geneID, len(coding_sequence)
	
			if DNA_distance >= 0:
				if options.mask:
					protein_distance = -1
				else:
					protein_distance = str_distance(str(pep), str(pepc))
	
				gene_model_analysis_file.write('\t' + str(DNA_distance))
				gene_model_analysis_file.write('\t' + str(len(str(pep))))
	
				if protein_distance >= 0:
					gene_model_analysis_file.write('\t' + str(protein_distance))
					gene_model_analysis_file.write('\t' + 'NA')
					gene_model_analysis_file.write('\t' + 'No')
				else:
					gene_model_analysis_file.write('\t' + 'NA')
					gene_model_analysis_file.write('\t' + 'NA')
					gene_model_analysis_file.write('\t' + 'Yes')
			else:
				gene_model_analysis_file.write('\t' + 'NA')
				gene_model_analysis_file.write('\t' + str(len(str(pep))))
				gene_model_analysis_file.write('\t' + 'NA')
				gene_model_analysis_file.write('\t' + str(len(str(coding_sequence_converted)) - len(str(coding_sequence))))
				gene_model_analysis_file.write('\t' + 'Yes')
	
			if len(options.annotation) > 0:
				if string.replace(geneID, '.v3.1', '') in gene_annotation.keys():
					gene_model_analysis_file.write('\t' + gene_annotation[string.replace(geneID, '.v3.1', '')])
	
				for index in range(len(experimental_datasets)):
					gene_model_analysis_file.write('\t' + dataset_gene_expression[experimental_datasets[index]][geneID])
	
			gene_model_analysis_file.write('\n')


#ORF_file.close()
ORF_converted_file.close()
gene_model_analysis_file.close()
intron_junction_file.close()

# preprocess data for visualization
intron_junction_file = open(prefix + '_intron_junction_sequence.txt', 'r')

intron_junction_count = {}

for line in intron_junction_file.readlines():
	sline = string.split(line)

	if sline[2] not in intron_junction_count.keys():
		intron_junction_count[sline[2]] = 0

	intron_junction_count[sline[2]] += 1

intron_junction_file.close()

intron_junction_counts_file = open(prefix + '_intron_junction_counts.txt', 'w')

intron_junction_counts_file.write('index' + '\t' + 'junction' + '\t' + 'counts' + '\n')

index = 1

for junction in intron_junction_count.keys():
	intron_junction_counts_file.write(str(index) + '\t' + junction + '\t' + str(intron_junction_count[junction]) + '\n')
	index += 1

intron_junction_counts_file.close()

# data visualization
print 'Data Visualization'
summary_file.write(time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + ' ' + 'Data Visualization' + '\n')

data_visualization_file = open(prefix + '.R', 'w')

data_visualization_file.write('library(ggplot2)' + '\n')
data_visualization_file.write('library(scales)' + '\n')
data_visualization_file.write('\n')
data_visualization_file.write('SNP = read.table(file="' + prefix + '_SNP_frequency.txt", header=T)' + '\n')
data_visualization_file.write('SNP = data.frame(SNP)' + '\n')
data_visualization_file.write('\n')
data_visualization_file.write('png(file="' + prefix + '_SNP_frequency.png", height=600, width=600)' + '\n')
data_visualization_file.write('ggplot(SNP, aes(frequency)) + geom_histogram(binwidth = 1) + scale_x_continuous(limits = c(0, 100))' + '\n')
data_visualization_file.write('dev.off()' + '\n')
data_visualization_file.write('\n')
data_visualization_file.write('indel = read.table(file="' + prefix + '_indel_frequency.txt", header=T)' + '\n')
data_visualization_file.write('indel = data.frame(indel)' + '\n')
data_visualization_file.write('\n')
data_visualization_file.write('png(file="' + prefix + '_indel_frequency.png", height=600, width=600)' + '\n')
data_visualization_file.write('ggplot(indel, aes(frequency)) + geom_histogram(binwidth = 1) + scale_x_continuous(limits = c(0, 100))' + '\n')
data_visualization_file.write('dev.off()' + '\n')
data_visualization_file.write('' + '\n')
data_visualization_file.write('CGA = read.table(file="' + prefix + '_candidate_gene_analysis.txt", header=T)' + '\n')
data_visualization_file.write('CGA = data.frame(CGA)' + '\n')
data_visualization_file.write('' + '\n')
data_visualization_file.write('png(file="' + prefix + '_coverage_histogram.png", height=600, width=600)' + '\n')
data_visualization_file.write('ggplot(CGA, aes(Coverage)) + geom_histogram(binwidth = 0.05)' + '\n')
data_visualization_file.write('dev.off()' + '\n')
data_visualization_file.write('' + '\n')
data_visualization_file.write('png(file="' + prefix + '_total_coverage_piechart.png", height=600, width=600)' + '\n')
data_visualization_file.write('pie.data <- data.frame(group = c("Coverage met", "Coverage not met"), value = c(' + str(total_coverage[1]) + ', ' + str(total_coverage[0]) + '))' + '\n')
data_visualization_file.write('blank_theme <- theme_minimal() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), panel.grid=element_blank(), axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold"))' + '\n')
data_visualization_file.write('ggplot(pie.data, aes(x="", y=value, fill=group)) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) + scale_fill_manual(name="", values=c("#E69F00", "#56B4E9")) + blank_theme + theme(axis.text.x=element_blank()) + geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)])), label = percent(c(' + str(float(total_coverage[1]) / sum(total_coverage)) + ',' + str(float(total_coverage[0]) / sum(total_coverage)) + ')), size=5)' + '\n')
data_visualization_file.write('dev.off()' + '\n')
data_visualization_file.write('' + '\n')

if len(intron_junction_count.keys()) > 0:
	data_visualization_file.write('intron = read.table(file="' + prefix + '_intron_junction_counts.txt", header=T)' + '\n')
	data_visualization_file.write('intron = data.frame(intron)' + '\n')
	data_visualization_file.write('\n')
	data_visualization_file.write('png(file="' + prefix + '_intron_junction_frequency.png", height=500, width=1500)' + '\n')
	data_visualization_file.write('ggplot(intron, aes(index, counts)) + geom_bar(stat="identity") + scale_x_discrete(breaks=intron$index, labels=intron$junction)' + '\n')
	data_visualization_file.write('dev.off()' + '\n')
	data_visualization_file.write('\n')

data_visualization_file.close()

commands.getstatusoutput('R --vanilla < ' + prefix + '.R')

summary_file.write('Run stop time: ' + time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + '\n')
summary_file.close()
