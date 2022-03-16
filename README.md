# Analysis of the human reads in gut metagenome shotgun sequencing data
This is a repository of the codes used in the Tomofuji et al (Reconstruction of the personal information from the contaminated human reads in the gut metagenome shotgun sequencing data).  
Our script recovers following information from the human reads in the metagenome shotgun sequencing data  
・Genetic sex of the metagenome shotgun sequencing data  
・Pair of the metagenome shotgun sequencing data and genotype data derived from the same individual  
・Genetic ancestry of the metagenome shotgun sequencing data   

# Overview
<div align="center">
<img src="Figure/Graphical_abstract.jpg" width=60%>
</div>

# Requirements
・R (version 4.0.1)  
・Trimmomatic (version 0.39)  
・tidyverse (version 1.3.0)  
・python3 (version 3.7.6)  
・Picard (version 2.22.8)  
・bowtie2 (version 2.3.5.1)     
・samtools (version 2.3.5.1)   
・bedtools (version 2.29.2)  
・bcftools (version 1.10.2)  
・fastqc (version 0.11.9)  
・plink (version 1.90b4.4)  

# 1. Extraction of the human reads and prediction of genetic sex
<div align="center">
<img src="Figure/human_read_extraction_figure.jpg" width=40%>
</div>

First, human reads are extracted from gut metagenome shotgun sequencing data with the script `PIPELINE_1_human_read_extraction.sh`.  
Input file should be named as `${ID}_R1.fastq.gz` and `${ID}_R2.fastq.gz`  
Following variables are required:  
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
<div align="center">
<img src="Figure/likelihood_score.jpg" width=100%>
</div>

Likelihood scores for each metagenome shotgun sequencing data are calculated with the script `PIPELINE_2_likelihood_reidentification_test.sh`. The likelihood score reflects the likelihood that the observed human reads in the gut MSS data are derived from the target genotype data in `GENOME_FILE`.   
Following variables are required:  
`VCF_NAME`: Name of the vcf file added to the output filename  
`GENOME_FILE`: Vcf file for genotype dataset for which the likelihood score is calculated  
`REF_FASTA`: reference genome file (ex: hg37_1kg_decoy)  
`DIR`: Directory for analysis  
`REFERENCE_PANEL`: Vcf file for reference genotype data  
`REF_PLINK`: Plink file for reference genotype data  
`ID`: Sample ID   

This script outputs likelihood score for each pair of the metagenome shotgun sequencing data and target genotype data.
The columns of the `${ID}_likelihood_p_val_result.txt` indicate following values

`ID`: Sample ID of the target genotype data   
`Score`: Lilelihood score  
`EMP_P`: Empilically caluclated P-values   
`ANA_p`: P-values analytically calculated from the standardized likelihood score  
`RANK`: Rank of the likelihood score among the genotype dataset  
`ID_MATCH`: Whether the ID of the metagenome shotgun sequencing data is matched to that of the target genotype data  

# 3. Prediction of the ancestry
<div align="center">
<img src="Figure/likelihood_score_anc.jpg" width=100%>
</div>

Likelihood scores for each metagenome shotgun sequencing data are calculated with the script `PIPELINE_3_likelihood_test_ancestry.sh`. The likelihood score reflects the likelihood that the observed human reads in the gut MSS data are derived from the specified ancestries (`AMR, AFR, EUR, EAS, SAS` in this study).  
Following variables are required:  

`REF_FASTA`: reference genome file (ex: hg37_1kg_decoy)   
`DIR`: Directory for analysis  
`OKG_REF_DIR`: Directory for the 1KG reference data   
`ID`: Sample ID   

Following files should be in `OKD_REF_DIR`   
1. ALL_POP.chr{1-22}.freq.chr.pos.gz

|  CHR  |  POS |  REF  |  ALT |  ID  |  AMR |  AFR  |  EUR |  EAS  |  SAS  |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|  22  |  16050075 |  A  |  G |  rs587697622  |  0 |  0  |  0 |  0  |  0.001022  |
|  22  |  16050115 |  G  |  A |  rs587755077  |  0.001441 |  0  |  0.02345 |  0  |  0 |
|  ...  |  ... |  ...  |  ... |  ...  |  ... |  ...  |  ... |  ...  |  ... |

`CHR`: Chromosome
`POS`: Position
`REF`: Reference allele
`ALT`: Alternative allele  
`ID`: SNP ID
`AMR`-`SAS`: Allele frequency in each ancestry


2. 1KG_CHR_POS_REF_ALT_sorted.txt (sorted based on the CHR:POS)

|  CHR:POS |  REF  |  ALT |
| --- | --- | --- |
|  chr10:100000003  |  C  |  T  |
|  chr10:10000001  |  C  |  T  |
|  ...  |  ... |  ...  | 

`CHR:POS`: Chromosome:Position
`REF`: Reference allele
`ALT`: Alternative allele  

The output is written to `${ID}_population_likelihood_result.txt` and each columun indicates following values

`Sample_ID`: Sample ID of the metagenome shotgun sequencing data  
`USED_BASE_NUM`: Number of the bases used for the analysis   
`{AMR, EUR, AFR, EAS, SAS}_LIK`: Likelihood score for each ancestry   
`TOP_POP`: The population with the highest likelihood score    


## Contact
Yoshihiko Tomofuji: ytomofuji_at_sg.med.osaka-u.ac.jp
