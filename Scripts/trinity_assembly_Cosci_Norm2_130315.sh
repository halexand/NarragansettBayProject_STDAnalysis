#!/bin/bash

#PBS -k o
#PBS -q shared
#PBS -l nodes=1:ppn=32,vmem=50gb,walltime=50:00:00
#PBS -M alexandh@mason.indiana.edu
#PBS -m abe
#PBS -N Trinity
#PBS -j oe

module load bowtie
module load java
module load trinityrnaseq/2013-02-25

ulimit -s unlimited

wd=/N/dc2/projects/ldeo/rawData/Coscinodiscus_AmazonIsolate/trimmed/normalized/

cd $wd


outputDir=Trinity_output/ #Added output directory for Trinity results
mkdir -p $outputDir

/N/soft/mason/trinityrnaseq/2013-02-25/Trinity.pl --seqType fq --JM 10G --left SH221_ATCACG_L005__trimmed_paired_norm_1.fastq SH222_CGATGT_L005__trimmed_paired_norm_1.fastq SH223_TTAGGC_L005__trimmed_paired_norm_1.fastq SH230_GTGGCC_L0017cat__trimmed_paired_norm_1.fastq --right SH221_ATCACG_L005__trimmed_paired_norm_2.fastq SH222_CGATGT_L005__trimmed_paired_norm_2.fastq SH223_TTAGGC_L005__trimmed_paired_norm_2.fastq SH230_GTGGCC_L0017cat__trimmed_paired_norm_2.fastq --CPU 10 --output $outputDir


exit 0
(END)
