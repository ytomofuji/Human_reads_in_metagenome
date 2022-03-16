#!/bin/sh
<<VARIABLE
VCF_NAME: Name of the vcf file added to the output files
GENOME_FILE: Vcf file for genotype dataset for which the likelihood score is calculated
REF_FASTA: reference genome file (ex: hg37_1kg_decoy)
DIR: Directory for analysis
REFERENCE_PANEL: Vcf file for reference genotype data
REF_PLINK: Plink file for reference genotype data
ID: Sample ID
VARIABLE

mkdir -p ${DIR}/LIKELIHOOD_PVAL_${VCF_NAME}/${ID}
cd ${DIR}/LIKELIHOOD_PVAL_${VCF_NAME}/${ID}

bcftools mpileup \
-f ${REF_FASTA} \
-Oz \
-I \
-q 40 \
-Q 20 \
-a FORMAT/AD,FORMAT/DP,FORMAT/DV,FORMAT/DP4,FORMAT/DPR \
-T ${GENOME_FILE} \
${DIR}/HUMAN_READ/${ID}/BAM/non_bacterial_${ID}_mapped.bam > non_bacterial_${ID}_mapped.pileup

bcftools view non_bacterial_${ID}_mapped.pileup | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP]\n' > \
non_bacterial_${ID}_mapped.depth


#Clumping based on the reference genome panel
bcftools view \
-T non_bacterial_${ID}_mapped.depth \
-v snps \
${REFERENCE_PANEL} | \
bcftools query -f '%ID\tChr%CHROM:%POS\n' \
> non_bacterial_${ID}_refpanel_SNP_list.txt

awk -F"\t" '{print $1}' non_bacterial_${ID}_refpanel_SNP_list.txt \
> non_bacterial_${ID}_refpanel_SNP_before_pruning


plink \
--threads 1 \
--memory 16000 \
--bfile ${REF_PLINK} \
--snps-only \
--autosome \
--indep-pairwise 100 30 0.1 \
--keep-allele-order \
--extract non_bacterial_${ID}_refpanel_SNP_before_pruning \
--out ${ID}.pruned.snp

rm non_bacterial_${ID}_refpanel_SNP_before_pruning

#calculate MAF
plink \
--threads 1 \
--memory 16000 \
--bfile ${REF_PLINK} \
--snps-only \
--autosome \
--keep-allele-order \
--extract ${ID}.pruned.snp.prune.in \
--freq \
--out ${ID}.pruned.snp

cat non_bacterial_${ID}_refpanel_SNP_list.txt | sort -k 1,1 \
> non_bacterial_${ID}_refpanel_SNP_list.sorted.txt

cat ${ID}.pruned.snp.frq | awk 'NR > 1 {print $0}' | sort -k 2,2 | \
join -1 1 -2 2 non_bacterial_${ID}_refpanel_SNP_list.sorted.txt - \
> ${ID}.pruned.snp.frq.merged
gzip -f ${ID}.pruned.snp.frq.merged
rm ${ID}.pruned.snp.frq

zcat ${ID}.pruned.snp.frq.merged.gz | awk '{print $2}' | sed 's/Chr//g' | sed 's/:/\t/g' \
> non_bacterial_${ID}_refpanel_SNP_CHR_POS

#extract from genotype dataset
bcftools view \
-R non_bacterial_${ID}_refpanel_SNP_CHR_POS \
-v snps \
${GENOME_FILE} | \
bcftools query -f 'Chr%CHROM:%POS\t%REF\t%ALT\n' | \
sort -k 1,1 \
> non_bacterial_${ID}_ref_allele_list.txt

awk -F'\t' '{print "Chr"$1":"$2"\t"$0}' non_bacterial_${ID}_mapped.depth | \
sort -k 1,1 -t$'\t' | join -j 1 -t$'\t' non_bacterial_${ID}_ref_allele_list.txt - \
> non_bacterial_${ID}_mapped.depth.merged

python3 \
/path/to/ALLELE_COUNT.py \
non_bacterial_${ID}_mapped.depth.merged \
non_bacterial_${ID}_mapped.alelle.count

bcftools view \
-R non_bacterial_${ID}_refpanel_SNP_CHR_POS \
-v snps \
${GENOME_FILE} | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed 'sX|X/Xg' | sed 'sX1/0X0/1Xg' \
> GENOTYPE.snp.truth.all.sample.bayes.${ID}

SAMPLE_ID=`bcftools query --list-samples \
${GENOME_FILE} | \
tr '\n' '\t' | \
sed 's/^/CHROM\tPOS\tREF\tALT\t/g' | sed 's/\t$//g'`

sed -i "1i ${SAMPLE_ID}" GENOTYPE.snp.truth.all.sample.bayes.${ID}
gzip -f GENOTYPE.snp.truth.all.sample.bayes.${ID}

# Calculate likelihood score
Rscript /path/to/RCODE_LIKELIHOOD_CALC_PVAL.R \
non_bacterial_${ID}_mapped.alelle.count \
GENOTYPE.snp.truth.all.sample.bayes.${ID}.gz \
${ID} ${ID} \
${ID}.pruned.snp.frq.merged.gz \
1e-6 \
1e_6 \
99999

gzip -f ${ID}_permutation_result.txt
gzip -f ${ID}_likelihood_p_val_summary.txt
gzip -f ${ID}_likelihood_p_val_result.txt

gzip -f non_bacterial_${ID}_mapped.alelle.count
gzip -f non_bacterial_${ID}_refpanel_SNP_CHR_POS
gzip -f ${ID}.pruned.snp.prune.in
gzip -f non_bacterial_${ID}_mapped.depth.merged

rm *.nosex
rm *.log
rm non_bacterial_${ID}_refpanel_SNP_list.txt
rm non_bacterial_${ID}_refpanel_SNP_list.sorted.txt
rm non_bacterial_${ID}_mapped.pileup
rm non_bacterial_${ID}_mapped.depth
rm .pruned.snp.prune.out
rm non_bacterial_${ID}_ref_allele_list.txt