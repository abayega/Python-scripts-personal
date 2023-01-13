#!/bin/bash

export PATH=/gs/project/wst-164-ab/anthony/software/bin:$PATH

cd ~/../../gs/project/wst-164-ab/anthony/Nanopore_data/C004_06_5_290916/nanook/fasta/pass/2D


cat pass/*.fasta | awk '//{print }' >combinedPass_C006_02_4.fasta

cat *.fasta | awk '//{print }' >combined2DPass_C002_05_4.fasta

awk 'BEGIN{y=0} !/^>/{y=y + length($0)} END{print y}' combined_PF_2D_C002_.fasta1

cat *.fasta | awk 'BEGIN{y=0} !/^>/{y=y + length($0); if ($0 >= 50000) print} END{print y}' > longRead_combined_PF_2D_C002_.fasta1

cat *.fasta | head -2 | awk '{if ($1 >= 50000) print} END{print $0:print $1}' > longRead_combined_PF_2D_C002_.fasta1

awk 
{
    n=index($0, ">")
    if (n > 0)
	n = NR
	#m = $0		
	y=0
	print($0)
    if ($1 ^[A-Z])
	y += length($1)
	x = NR
    if (n > x)
	print(m,y)


sed 's/5-1582/Scaffold4/g' Contig5-1582_bwa-mem.sam > Contig5-1582-Scassfold4_bwa-mem.sam
awk ' NR == 2418 {print(">Scaffold4\n"$0)}' ../connected_contigs.fasta > 2.txt
awk 'NR<6 && NR%2==0 {print(length($0))}' Contig5-1582.fasta

	

	

