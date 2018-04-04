#! /usr/bin/env python
# -*- coding: UTF-8 -*-

"""
QKgenome_missingdata.py [options] datasets coverage prefix
Author: Matthew Moscou <matthew.moscou@tsl.ac.uk>
Remove sites with extensive missing data in a Phylip formatted alignment of polymorphic SNP sites based on output from QKgenome_conversion.py
This performs the following:
	1. import Phylip formatted file
	2. export Phylip formatted file with removal of sites with extensive missing data

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
parser.add_option("-m", "--missing", action="store", dest="missing", default=100.0, help="Percent of missing data tolerated (between 0 and 100)")
parser.add_option("-a", "--alignment", action="store", dest="alignment", default="", help="Alignment file in Phylip format")
parser.add_option("-o", "--output", action="store", dest="output", default="", help="Output alignment in Phylip format after missing data removal")
(options, args) = parser.parse_args()

# read raw alignment
phylip_tree = open(options.alignment, 'r')

ID_alignment = {}

truth = False

for line in phylip_tree.readlines():
	sline = string.split(line)

	if truth:
		ID_alignment[sline[0]] = sline[1]

	truth = True

phylip_tree.close()

# generate alignment using user-specified threshold
ID_threshold = {}

for ID in ID_alignment.keys():
	ID_threshold[ID] = ''

for index in range(len(ID_alignment[ID_alignment.keys()[0]])):
	site_calls = []

	for ID in ID_alignment.keys():
		site_calls.append(ID_alignment[ID][index])

	missing_data = ((site_calls.count('N') * 1.0) / len(site_calls)) * 100.0

	if missing_data < float(options.missing):
		for ID in ID_alignment.keys():
			ID_threshold[ID] += ID_alignment[ID][index]


# export multiple sequence alignment of polymorphic sites
print 'export multiple sequence alignment of polymorphic sites'
phylip_input = open(options.output, 'w')

phylip_input.write(' ' + str(len(ID_threshold.keys())) + ' ' + str(len(ID_threshold[ID_threshold.keys()[0]])) + '\n')

for ID in ID_threshold.keys():
	phylip_input.write(ID + ''.join([' '] * (10 - len(ID))) + ID_threshold[ID] + '\n')

phylip_input.close()
