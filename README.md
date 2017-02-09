# QKgenome
A set of scripts for converting genomes based on resequencing information.

## Scripts
<i>QKgenome_conversion.py</i> takes as input a genome/region/gene sequence(s), GFF3, annotation, expression, and VarScan output for pileup2snp and pileup2indel and bedtools genomecov, processes the information based on user specified parameters for coverage and variant frequency, then exports a converted genome (relative to the reference) along with coverted sequence, GFF3, and CDS. 

<i>QKgenome_phylogeny.py</i> uses the output generated from <i>QKgenome_conversion.py</i> to generate an input file in Phylip format for generating a phylogenetic tree. In addition, the data can be used for association genetics.

## Example
Coming soon.
