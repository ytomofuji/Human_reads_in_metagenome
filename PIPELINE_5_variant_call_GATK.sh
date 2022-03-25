#!/bin/sh
<<VARIABLE
REF_FASTA: reference genome file (ex: hg37_1kg_decoy)
DIR: Directory for analysis
DBSNP_FILE: dbsnp file 
ID: Sample ID
VARIABLE

mkdir -p ${DIR}/GATK/${ID}
cd ${DIR}/GATK/${ID}

#Variant call with GATK
/path/to/GATK4/gatk \
--java-options "-Xmx12g" \
HaplotypeCaller \
-R ${REF_FASTA} \
-I ${DIR}/HUMAN_READ/${ID}/BAM/non_bacterial_${ID}_mapped.bam \
--dbsnp ${DBSNP_FILE} \
-O Unfiltered_GATK_${ID}.vcf

#make snp only data set
for i in `seq 1 22`
do
echo -e ${i} >> auto.list
done

#only autosome
/path/to/GATK4/gatk \
--java-options "-Xmx12g" \
SelectVariants \
-V Unfiltered_GATK_${ID}.vcf \
-select-type SNP \
-L auto.list \
-O Auto_Unfiltered_GATK_${ID}.vcf

rm auto.list
rm Unfiltered_GATK_${ID}.vcf

#hard filtering on snp dataset
/path/to/GATK4/gatk \
--java-options "-Xmx8g" \
VariantFiltration \
-V Auto_Unfiltered_GATK_${ID}.vcf \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "DP < 2.0" --filter-name "DP2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O GATK_${ID}.filtered.vcf

bcftools view -i "%FILTER='PASS'" GATK_${ID}.filtered.vcf \
> GATK_${ID}.vcf
rm GATK_${ID}.filtered.vcf