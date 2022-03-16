# Analysis of the human reads in gut metagenome shotgun sequencing data
This is a repository of the codes used in the Tomofuji et al (Reconstruction of the personal information from the contaminated human reads in the gut metagenome shotgun sequencing data).
Our software recovers following information from the human reads in the metagenome shotgun sequencing data
・Genetic sex of the metagenome shotgun sequencing data
・Pair of the metagenome shotgun sequencing data and genotype data derived from the same individual
・Genetic ancestry of the metagenome shotgun sequencing data

# Overview
![Image 1](https://github.com/ytomofuji/Human_reads_in_metagenome/blob/main/Figure/Graphical_abstract.jpg)

# Requirements
```
bowtie2  
samtools  
bedtools  
bcftools  
fastqc  
```

# 1. Extraction of the human reads + prediction of genetic sex
First, human reads are extracted from gut metagenome shotgun sequencing data with the script `PIPELINE_1_human_read_extraction.sh`. 
Following variables are required

`TRIM_ADAPT`: adaptor-sequence (ex: TruSeq3-PE-2.fa)  
`BOWTIE2_REF_HUM`: reference file for bowtie2 (human, ex: hg37_1kg_decoy)  
`BOWTIE2_REF_BAC`: reference file for bowtie2 (bacteria)  
`DIR`: Directory for analysis  
`FASTQ_DIR`: Directory of original fastq file  
`ID`: Sample ID   
`BED_OF_NONPAR`: bed file for non-pseudoautosomal region (non-PAR) of the X and Y chromosomes  

This script outputs human reads in the metagenome shotgun sequencing data (`non_bacterial_${ID}_mapped.bam`) which can be used in the subsequent analyses.
In addition, coverages of the non-PAR of X and Y chromosomes are output into `bedcov_${ID}_nonPAR_XY.txt`, which can be used to predict genetic sex of the metagenome shotgun sequencing data.

# 2. Re-identification from a set of genotype data
Pair of the metagenome shotgun sequencing data and genotype data derived from the same individual
