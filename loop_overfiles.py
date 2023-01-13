#!/bin/env python3

import os

def loop_files(dir,pass_or_fail):
	k = pass_or_fail
	files = os.listdir(dir)
	for file in files:
		if file.startswith('201709'):
			filename = file.split('_')
			filename = '_'.join([filename[2],filename[3],filename[4],filename[1]])
			dirname = dir+'/'+file			
			os.system('cat %s/workspace/%s/*.fastq > %s/workspace/%s/%s_%s.fq' %(dirname,k,dirname,k,filename,k))
			os.system('rm %s/workspace/%s/*.fastq' %(dirname,k))
			
def loop_files2(dir,pass_or_fail):
	k = pass_or_fail
	files = os.listdir(dir)
	for file in files:
		if file.isdecimal():
			filename = file.split('_')
			#filename = '_'.join([filename[2],filename[3],filename[4],filename[1]])
			filename = 'C010_04_4_1350'
			dirname = dir+'/'+file			
			os.system('cat %s/workspace/%s/*.fastq > %s/workspace/%s/%s_%s.fq' %(dirname,k,dirname,k,filename,k))
			os.system('rm %s/workspace/%s/*.fastq' %(dirname,k))
			
#loop_files('/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/reads/nanopore/capitata/rna_seq','fail')
loop_files2('/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/reads/nanopore/capitata/rna_seq/20170916_1350_C010_04_4/basecalled','pass')
#loop_files('/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/reads/nanopore/capitata/rna_seq/pass','fail')