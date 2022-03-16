#!/bin/sh
<<VARIABLE
REF_FASTA: reference genome file (ex: hg37_1kg_decoy)
DIR: Directory for analysis
OKG_REF_DIR: Directory of the 1KG files
ID: Sample ID
VARIABLE

mkdir -p ${DIR}/LIKELIHOOD_POP_1KG/${ID}
cd ${DIR}/LIKELIHOOD_POP_1KG/${ID}

bcftools mpileup \
-f ${REF_FASTA} \
-Oz \
-I \
-q 40 \
-Q 20 \
-a FORMAT/AD,FORMAT/DP,FORMAT/DV,FORMAT/DP4,FORMAT/DPR \
-T ${OKG_REF_DIR}/1KG_CHR_POS_list.txt \
${DIR}/HUMAN_READ/${ID}/BAM/non_bacterial_${ID}_mapped.bam > non_bacterial_${ID}_mapped.pileup

bcftools view non_bacterial_${ID}_mapped.pileup | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD\t%DP]\n' > \
non_bacterial_${ID}_mapped.depth

awk -F'\t' '{print "Chr"$1":"$2"\t"$0}' non_bacterial_${ID}_mapped.depth | \
sort -k 1,1 -t$'\t' | join -j 1 -t$'\t' ${OKG_REF_DIR}/1KG_CHR_POS_REF_ALT_sorted.txt - \
> non_bacterial_${ID}_mapped.depth.merged

python3 \
/path/to/ALLELE_COUNT.py \
non_bacterial_${ID}_mapped.depth.merged \
non_bacterial_${ID}_mapped.alelle.count

gzip -f non_bacterial_${ID}_mapped.alelle.count

Rscript \
/path/to/RCODE_LIKELIHOOD_POP_CALC.R \
non_bacterial_${ID}_mapped.alelle.count.gz \
${OKG_REF_DIR} \
${ID}

rm non_bacterial_${ID}_mapped.depth
rm non_bacterial_${ID}_mapped.pileup
rm non_bacterial_${ID}_mapped.depth.merged
