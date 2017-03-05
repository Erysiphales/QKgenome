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

### Optional commands
1. Annotation

   -a, --annotation   To include a tab-delimited file containing gene annotations
2. Masking of regions without read coverage

   -m, --mask         Will convert nucleotides to N where read coverage is below the provided threshold
3. Inclusion of heterozygous SNPs

   Requires two parameters:
   
      -l, --lower     Lower bound for SNP frequency threshold
      
      -u, --upper     Upper bound for SNP frequency threshold

<i>QKgenome_phylogeny.py</i> uses the output generated from <i>QKgenome_conversion.py</i> to generate an input file in Phylip format for generating a phylogenetic tree. Additional uses of data include association genetics or candidate gene analysis.

Note: Some customization is required to take full advantage of <i>QKgenome_phylogeny.py</i>. In particular, it is a common feature of genomes that multiple gene models exist for a gene. These arise from splice site variants, which can create transcripts with different 5' and 3' UTRs, as well as change the coding sequence. There is not a single accepted standard nomenclature for gene models, therefore, removal of potential redundancy in an analysis requires modification of the script for the species specific nomenclature. There is existing code within the script to handle this, but the user needs to modify it.

## Software and modules
The following software and modules are requires to run the suite of QKgenome scripts:

1. Python 2.7
  * BioPython 1.64
2. R 3.2.3
  * ggplots2
  * scales

Other versions are likely to be functional, but the versions described above were used in the development of the scripts.

## Example
### Convert a genomic region in <i>Brachypodium distachyon</i> from reference sequence based on a resequenced genotype
Converts a region of the <i>Brachypodium distachyon</i> genome from the reference (accession Bd21) into the accession ABR6 version. Read coverage required is 20 reads and a minimum variant frequency of reads at 95%. In addition, read counts from RNAseq data are included within the analysis.
```bash
python QKgenome_conversion.py 20 95.0 Bd4_24522081_29856080.fa Bd4_24522081_29856080.gff3 Yrr1_Jer1_sorted.rmdup.pileup2snp.txt Yrr1_Jer1_sorted.rmdup.pileup2indel.txt Yrr1_Jer1_sorted.rmdup.genomecov.txt Yrr1_Jer1 Yrr1_Jer1_reference_Jer1_RNAseq_tophat_readCounts.txt
```

### Convert the <i>Puccinia striiformis</i> f. sp. <i>tritici</i> genome based on resequenced isolates from the US, UK, and Australia
```bash
python QKgenome_conversion.py 20 90.0 puccinia_striiformis_pst-78_1_supercontigs.fasta puccinia_striiformis_pst-78_1_transcripts.gff3 PST_PST_104E137A-_RNAseq_sorted.rmdup.pileup2snp.txt PST_PST_104E137A-_RNAseq_sorted.rmdup.pileup2indel.txt PST_PST_104E137A-_RNAseq_sorted.genomecov.txt PST_78_PST_104E137A-_RNAseq

python QKgenome_conversion.py 20 90.0 puccinia_striiformis_pst-78_1_supercontigs.fasta puccinia_striiformis_pst-78_1_transcripts.gff3 PSH_B012_RNAseq_trimmomatic_sorted.pileup2snp.txt PSH_B012_RNAseq_trimmomatic_sorted.pileup2indel.txt PSH_B012_RNAseq_trimmomatic_sorted.genomecov.txt PST_78_PSH_B012_RNAseq

python QKgenome_conversion.py 20 90.0 puccinia_striiformis_pst-78_1_supercontigs.fasta puccinia_striiformis_pst-78_1_transcripts.gff3 PST_0821_gDNA_trimmomatic_sorted.rmdup.pileup2snp.txt PST_0821_gDNA_trimmomatic_sorted.rmdup.pileup2indel.txt PST_0821_gDNA_trimmomatic_sorted.rmdup.genomecov.txt PST_78_PST_0821_gDNA

python QKgenome_conversion.py 20 90.0 puccinia_striiformis_pst-78_1_supercontigs.fasta puccinia_striiformis_pst-78_1_transcripts.gff3 PST_08501_gDNA_trimmomatic_sorted.rmdup.pileup2snp.txt PST_08501_gDNA_trimmomatic_sorted.rmdup.pileup2indel.txt PST_08501_gDNA_trimmomatic_sorted.rmdup.genomecov.txt PST_78_PST_08501_gDNA

python QKgenome_conversion.py 20 90.0 puccinia_striiformis_pst-78_1_supercontigs.fasta puccinia_striiformis_pst-78_1_transcripts.gff3 PST_1108_gDNA_trimmomatic_sorted.rmdup.pileup2snp.txt PST_1108_gDNA_trimmomatic_sorted.rmdup.pileup2indel.txt PST_1108_gDNA_trimmomatic_sorted.rmdup.genomecov.txt PST_78_PST_1108_gDNA

python QKgenome_conversion.py 20 90.0 puccinia_striiformis_pst-78_1_supercontigs.fasta puccinia_striiformis_pst-78_1_transcripts.gff3 PST_1108_RNAseq_trimmomatic_sorted.pileup2snp.txt PST_1108_RNAseq_trimmomatic_sorted.pileup2indel.txt PST_1108_RNAseq_trimmomatic_sorted.genomecov.txt PST_78_PST_1108_RNAseq

python QKgenome_conversion.py 20 90.0 puccinia_striiformis_pst-78_1_supercontigs.fasta puccinia_striiformis_pst-78_1_transcripts.gff3 PST_78_gDNA_trimmomatic_sorted.rmdup.pileup2snp.txt PST_78_gDNA_trimmomatic_sorted.rmdup.pileup2indel.txt PST_78_gDNA_trimmomatic_sorted.rmdup.genomecov.txt PST_78_PST_78_reference

python QKgenome_phylogeny.py PST_78_analysis_inventory.txt 1.0 PST_78_analysis_information.txt PST_78_analysis_1.0.phy
```
