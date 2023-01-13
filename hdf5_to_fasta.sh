#!/bin/bash

h5dump -d "/Analyses/Basecall_1D_000/BaseCalled_template/Fastq" *.fast5 | awk 'BEGIN{print ">"} NR==11,NR==12{print}' > file.fasta

/Analyses/Basecall_1D_000/BaseCalled_template/Events

h5dump -d "/Analyses/Basecall_1D_000/BaseCalled_template/Events" tiare_20160915_FNFAB29949_MN18767_mux_scan_C002_05_6_23523_ch12_read77_strand.fast5 | awk 'BEGIN{print ">"} //{print}' > file.fasta
