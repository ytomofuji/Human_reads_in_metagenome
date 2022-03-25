#!/bin/sh
<<VARIABLE
REF_FASTA: reference genome file (ex: hg37_1kg_decoy)
DIR: Directory for analysis
REF_PANEL_DIR: Directory of the 1KG files (per chromosome genome file is Chr{1-22}.vcf.gz)
MAP_FILE_DIR: Directory of the 1KG files (per chromosome map file is plink.chr${1-22}.GRCh37.map)
ID: Sample ID
VARIABLE

mkdir -p ${DIR}/ldWGS/${ID}
cd ${DIR}/ldWGS/${ID}

echo "${ID}" > ${ID}.new.sample.name

for CHR in `seq 1 22`
do
#variant call by samtools
bcftools mpileup --regions ${CHR} \
-f ${REF_FASTA} \
-I \
-q 20 \
-Q 20 \
-Ou ${DIR}/HUMAN_READ/${ID}/BAM/non_bacterial_${ID}_mapped.bam | \
bcftools call --skip-variants indels -c -Oz -v -o Unfiltered_Samtools_${ID}_${CHR}.vcf.gz

#2 step imputation
java -Xss5m -Xmx16g -jar \
/path/to/beagle4/beagle.27Jan18.7e1.4.1.jar \
gl=Unfiltered_Samtools_${ID}_${CHR}.vcf.gz \
ref=${REF_PANEL_DIR}/Chr${CHR}.vcf.gz \
out=${ID}_${CHR}.bgl.4.1 \
map=${MAP_FILE_DIR}/plink.chr${CHR}.GRCh37.map \
niterations=0 \
modelscale=2 \
nthreads=4

rm ${ID}_${CHR}.bgl.4.1.log

java -Xss5m -Xmx16g -jar \
/path/to/beagle5/beagle.18May20.d20.5.1.jar \
gt=${ID}_${CHR}.bgl.4.1.vcf.gz \
ref=${REF_PANEL_DIR}/Chr${CHR}.vcf.gz \
out=${ID}_${CHR}.bgl.5.1 \
map=${MAP_FILE_DIR}/plink.chr${CHR}.GRCh37.map \
nthreads=4 \
gp=true

rm ${ID}_${CHR}.bgl.5.1.log

rm Unfiltered_Samtools_${ID}_${CHR}.vcf.gz
rm ${ID}_${CHR}.bgl.4.1.vcf.gz

bcftools reheader -s ${ID}.new.sample.name \
${ID}_${CHR}.bgl.5.1.vcf.gz > ${ID}_${CHR}.bgl.5.1.reheadered.vcf.gz

rm ${ID}_${CHR}.bgl.5.1.vcf.gz
rm ${ID}_${CHR}.new.sample.name
bcftools index -f -t ${ID}_${CHR}.bgl.5.1.reheadered.vcf.gz
done

rm ${ID}.new.sample.name