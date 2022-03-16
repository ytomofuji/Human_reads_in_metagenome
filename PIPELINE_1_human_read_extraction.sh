#!/bin/sh
<<VARIABLE
TRIM_ADAPT: adaptor-sequence (ex: TruSeq3-PE-2.fa)
BOWTIE2_REF_HUM: reference file for bowtie2 (human, ex: hg37_1kg_decoy)
BOWTIE2_REF_BAC: reference file for bowtie2 (bacteria)
DIR: Directory for analysis
FASTQ_DIR: Directory of original fastq file
ID: Sample ID
BED_OF_NONPAR: bed file for non-PAR
VARIABLE

#Requirement
# bowtie2, samtools, bedtools, fastqc

mkdir -p ${DIR}/HUMAN_READ/${ID}/FASTQ
cd ${DIR}/HUMAN_READ/${ID}/FASTQ

#read QC
java -Xmx16g -jar \
/path/to/trimmomatic/trimmomatic.jar PE -threads 1 -phred33 \
${FASTQ_DIR}/${ID}_R1.fastq.gz \
${FASTQ_DIR}/${ID}_R2.fastq.gz \
paired_${ID}_R1.fastq.gz unpaired_${ID}_R1.fastq.gz \
paired_${ID}_R2.fastq.gz unpaired_${ID}_R2.fastq.gz \
ILLUMINACLIP:${TRIM_ADAPT}:2:30:10:8:true \
TRAILING:20 MINLEN:60

cd ${DIR}/HUMAN_READ/${ID}

mkdir -p BAM
cd BAM

NUM0=`zcat ${FASTQ_DIR}/${ID}_R1.fastq.gz | wc -l | awk '{print $1/4}'`
NUM1=`zcat ../FASTQ/paired_${ID}_R1.fastq.gz | wc -l | awk '{print $1/4}'`

bowtie2 -R 3 --no-discordant \
-x ${BOWTIE2_REF_HUM} \
-1 ../FASTQ/paired_${ID}_R1.fastq.gz \
-2 ../FASTQ/paired_${ID}_R2.fastq.gz \
-S ${ID}.sam

#extract only properly mapped read
samtools view -Sb -@ 2 ${ID}.sam | samtools sort -@ 2 -m 20G | \
samtools view -u -f 3 > ${ID}_mapped_pre.bam
rm ${ID}.sam
samtools index ${ID}_mapped_pre.bam

#add read group "@RG\tID:FLOWCELLID\tSM:${ID}\tPL:illumina\tLB:${ID}_library_1"
java -Xmx16g -jar /path/to/picard/picard.jar \
AddOrReplaceReadGroups I=${ID}_mapped_pre.bam \
O=${ID}_mapped.bam \
RGID=FLOWCELLID \
RGLB=${ID}_library_1 \
RGPU=DUMMY \
RGPL=illumina \
RGSM=${ID} \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=true

NUM2=`samtools view -c ${ID}_mapped.bam | awk '{print $1/2}'`

#remove duplicates by Picard
java -Xmx16g -jar /path/to/picard/picard.jar \
MarkDuplicates I=${ID}_mapped.bam \
O=rm_dup_${ID}_mapped.bam \
M=rm_dup_${ID}_mapped.metrics.txt \
REMOVE_DUPLICATES=true \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=LENIENT

samtools index rm_dup_${ID}_mapped.bam

NUM3=`samtools view -c rm_dup_${ID}_mapped.bam | awk '{print $1/2}'`

#extract only-non-bacterial reads
java -Xmx16g -jar /path/to/picard/picard.jar \
SamToFastq I=rm_dup_${ID}_mapped.bam \
F=rm_dup_${ID}_mapped_R1.fastq F2=rm_dup_${ID}_mapped_R2.fastq

gzip -f rm_dup_${ID}_mapped_R1.fastq
gzip -f rm_dup_${ID}_mapped_R2.fastq

bowtie2 -x \
${BOWTIE2_REF_BAC} \
-1 rm_dup_${ID}_mapped_R1.fastq.gz \
-2 rm_dup_${ID}_mapped_R2.fastq.gz \
-S non_bacterial_${ID}.sam

samtools view -b -f 12 -F 256 non_bacterial_${ID}.sam \
> non_bacterial_${ID}.bam

rm non_bacterial_${ID}.sam
samtools index non_bacterial_${ID}.bam

NUM4=`samtools view -c non_bacterial_${ID}.bam | awk '{print $1/2}'`

samtools view non_bacterial_${ID}.bam | \
cut -f 1 | sort -k 1,1 | uniq > non_bacterial_${ID}.txt

java -Xmx16g -jar /path/to/picard/picard.jar \
FilterSamReads I=rm_dup_${ID}_mapped.bam \
O=non_bacterial_${ID}_mapped.bam \
RLF=non_bacterial_${ID}.txt \
FILTER=includeReadList

samtools index non_bacterial_${ID}_mapped.bam

#retain only non_bacterial_${ID}_mapped.bam
rm non_bacterial_${ID}.bam
rm rm_dup_${ID}_mapped.bam
rm ${ID}_mapped_pre.bam
rm ${ID}_mapped.bam
rm non_bacterial_${ID}.txt
rm rm_dup_${ID}_mapped.bam.bai
rm non_bacterial_${ID}.bam.bai
rm ${ID}_mapped_pre.bam.bai
rm ${ID}_mapped.bai

cd ..

mkdir -p metrics
cd metrics

#count following num
#1. num original read
#2. num QCed read
#3. num properly mapped read
#4. num dum-removed read
#5. num decontaminated read
echo -e "Original\tRead_QC\tMapping\tRemove_dup\tRemove_BAC" > BT2_read_number.txt
echo -e "${NUM0}\t${NUM1}\t${NUM2}\t${NUM3}\t${NUM4}" >> BT2_read_number.txt

samtools idxstats ../BAM/non_bacterial_${ID}_mapped.bam \
> bami_non_bacterial_${ID}_unique_mapped.txt

samtools view -u \
../BAM/non_bacterial_${ID}_mapped.bam | \
bedtools coverage -a ${BED_OF_NONPAR} \
-b stdin -counts -F 0.1 | \
awk 'BEGIN{OFS="\t"}{print $0"\t"$4/($3-$2)"\t"$4}' | \
sed '1i CHR\tSTART\tEND\tBP\tDepth\tCount' > bedcov_${ID}_nonPAR_XY.txt

mkdir -p ${DIR}/HUMAN_READ/FASTQC
fastqc -f bam -t 2 --nogroup ${DIR}/HUMAN_READ/${ID}/BAM/non_bacterial_${ID}_mapped.bam -o ${DIR}/HUMAN_READ/FASTQC

rm -r ${DIR}/HUMAN_READ/${ID}/FASTQ
rm ${DIR}/HUMAN_READ/${ID}/BAM/rm_dup_${ID}_mapped_R1.fastq.gz
rm ${DIR}/HUMAN_READ/${ID}/BAM/rm_dup_${ID}_mapped_R2.fastq.gz