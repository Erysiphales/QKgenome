# QKgenome
A set of scripts for converting genomes based on resequencing information.

## Scripts
<i>QKgenome_conversion.py</i> converts a reference sequence based on read alignments from a different accession/ecotype/cultivar/isolate based on user specified parameters for coverage and variant frequency. The reference sequence can be a genome, a set of genes, or transcriptome. Memory requirements for the script are proportional to the length of the sequence space.

### Input
1. Sequence (FASTA)
2. GFF3
3. Variant calls (from `VarScan pileup2snp` and `VarScan pileup2indel`)
4. Read coverage (from `bedtools genomecov`)
5. Optional: Annotation (tab-delimited file)
6. Optional: Expression read counts (from `featureCounts`)

### Output
1. Converted sequence (FASTA)
2. Converted GFF3
3. Converted coding sequence
4. Summary of coverage, SNPs, and nonsynonymous mutations; optional: expression levels and annotation
5. Intron junctions in reference and converted sequence
6. Figures displaying read coverage on gene space, SNP and InDel variant frequency above coverage threshold, and the frequency of intron splice junctions.

## Optional commands
1. Annotation
   -a, --annotation
   To include a tab-delimited file containing gene annotations
2. Masking of regions without read coverage
   -m, --mask
   Will convert nucleotides to N where read coverage is below the provided threshold
3. Inclusion of heterozygous SNPs
   Requires two parameters:
      -l, --lower
      Lower bound for SNP frequency threshold
      -u, --uper
      Upper bound for SNP frequency threshold

<i>QKgenome_phylogeny.py</i> uses the output generated from <i>QKgenome_conversion.py</i> to generate an input file in Phylip format for generating a phylogenetic tree. Additional uses of data include association genetics or candidate gene analysis.

## Example
Converts a region of the <i>Brachypodium distachyon</i> genome from the reference (accession Bd21) into the accession ABR6 version. Read coverage required is 20 reads and a minimum variant frequency of reads at 95%. In addition, read counts from RNAseq data are included within the analysis.
```bash
python QKgenome_conversion.py 20 95.0 Bd4_24522081_29856080.fa Bd4_24522081_29856080.gff3 Yrr1_Jer1_sorted.rmdup.pileup2snp.txt Yrr1_Jer1_sorted.rmdup.pileup2indel.txt Yrr1_Jer1_sorted.rmdup.genomecov.txt Yrr1_Jer1 Yrr1_Jer1_reference_Jer1_RNAseq_tophat_readCounts.txt
```
