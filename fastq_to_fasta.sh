#! /bin/bash

gunzip -c file.fq.gz | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > file.fa

awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' read2.cutadapt.output.fastq > read2.cutadapt.output.fasta
#Get fastq read lengths
awk '{if(NR%4==2) {print(length($0))}}' Bo.E.5H_all_rnaseq_combined_fail.fastq > Bo.E.5H_all_rnaseq_combined_fail.readlengths
