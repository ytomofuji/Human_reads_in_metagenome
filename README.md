# Analysis of the human reads in gut metagenome shotgun sequencing data
This is a repository of the codes used in the Tomofuji et al (Reconstruction of the personal information from the contaminated human reads in the gut metagenome shotgun sequencing data).
Our software recovers following information from the human reads in the metagenome shotgun sequencing data
・Genetic sex of the metagenome shotgun sequencing data
・Pair of the metagenome shotgun sequencing data and genotype data derived from the same individual
・Genetic ancestry of the metagenome shotgun sequencing data

# Overview
![Image 1](https://github.com/ytomofuji/Human_reads_in_metagenome/blob/main/Graphical_abstract.jpg)

# Requirements
・bowtie2
・samtools
・bedtools
・bcftools
・fastqc

# 1. Extraction of the human reads
First, we extracted human reads from gut metagenome shotgun sequencing data with the script `PIPELINE_1_human_read_extraction.sh`. 
