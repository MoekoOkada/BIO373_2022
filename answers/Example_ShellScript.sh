#!/bin/sh
# HaplotypeCaller.sh
# Project: BIO373
# Date: 28-Sep-2022
# Author: Moeko Okada

## SOFTWARES
source /usr/local/ngseq/etc/lmod_profile
module load Aligner/BWA/0.7.17
module load Tools/samtools/1.11
module load Variants/GATK/4.2.0.0

bio373=/sratch/bio373_2022
workdir=${bio373}/moeko

cd ${workdir}
pwd

# indexing
bwa index 00_input/MedtrChr2.fa
samtools faidx 00_input/MedtrChr2.fa
gatk CreateSequenceDictionary -R 00_input.MedtrChr2.fa

# mapping
fMAP=0
if [ ${fMAP} -eq 1 ]; then
    for acc in "516950" "660389"; do
        bwa mem -M -t 2 -R "@RG\tID:CAV90ANXX.6\tPL:Illumina\tLB:${acc}\tSM:${acc}" \
            00_input/MedtrChr2.fa 00_input/${acc}_chr2_R{1,2}.fastq.gz | samtools sort -m 16G \
            -T /scratch/bio373_2021/moeko -o 01_aligned/${acc}.sorted.bam
    done
fi

# deduplication
fDEDUP=0
if [ ${fDEDUP} -eq 1 ]; then
    for acc in "516950" "660389"; do
        gatk MarkDuplicates -I 01_aligned/${acc}.sorted.bam \
            -O 02_dedup/${acc}.dedup.bam \
            -M 02_dedup/${acc}.metrics
        samtools index 02_dedup/${acc}.dedup.bam
    done
fi