# QKgenome
A set of scripts for converting genomes based on resequencing information.

## Scripts
<i>QKgenome_conversion.py</i> converts a reference sequence based on read alignments from a different accession/ecotype/cultivar/isolate based on user specified parameters for coverage and variant frequency. The reference sequence can be a genome, a set of genes, or transcriptome. Memory requirements for the script are proportional to the length of the sequence space.

### Input
1. Sequence (FASTA)
2. GFF3
3. Annotation (tab-delimited file)
4. Variant calls (from `VarScan pileup2snp` and `VarScan pileup2indel`)
5. Read coverage (from `bedtools genomecov`)
6. Optional: Expression read counts (from `featureCounts`)

### Output
1. Converted sequence (FASTA)
2. Converted GFF3
3. Converted coding sequence
4. Summary of coverage, SNPs, nonsynonymous mutations, and expression levels (optional)
5. Intron junctions in reference and converted sequence

<i>QKgenome_phylogeny.py</i> uses the output generated from <i>QKgenome_conversion.py</i> to generate an input file in Phylip format for generating a phylogenetic tree. Additional uses of data include association genetics or candidate gene analysis.

## Example

