#!/bin/bash

##Job settings
#PBS -V
#PBS -N HG002_alignment

##Job configuration
#PBS -l select=1:ncpus=64:mem=16gb
#PBS -lwalltime=10:00:00

#activate anaconda module
module load anaconda3/personal
source activate alignment

#set paths
DATADIR=/rds/general/user/fmazzaro/home/WORK/Large_variants/data/hg002_wes
REFGEN_OLD=/rds/general/project/lms-ware-raw/live/resources/reference/UCSC_hg19/allchrom.Chr1ToChrM.validated.fa
REFGEN_NEW=/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta

echo "Aligning HG002..."

fastqheader=$(zcat /rds/general/user/fmazzaro/home/WORK/Large_variants/data/hg002_wes/HG002.novaseq.wes_truseq.100x.R1.fastq.gz | head -n 1)
ID=$(echo $fastqheader | cut -d':' -f1,2,3,4 | cut -c2-)

bwa mem -Y -R @RG\\tPL:ILLUMINA\\tID:$ID\\tSM:$SAMPLEID -t 16 $REFGEN_NEW $DATADIR/HG002.novaseq.wes_truseq.100x.R1.fastq.gz $DATADIR/HG002.novaseq.wes_truseq.100x.R2.fastq.gz | samtools sort --threads 4 -o $DATADIR/HG002.bam

samtools index -@8 $DATADIR/HG002.bam