#!/usr/bin/env python

import re, os, sys

def rename_fasta_seq(fasta):
	ofasta = open(fasta)
	rfasta = ofasta.readlines()
	newfile = open('/gs/project/wst-164-ab/anthony/Nanopore_data/olive_fly/TULIP_assembly/olivefly_illumina_seeds_100pb.fasta', 'w')
	n = 1
	l = len(rfasta)
	for i in range(len(rfasta)):
		m = i + 1
		newfile.write('>%s\n' % m)
		newfile.write(rfasta[n])
		n += 2
		if n >= l:
			break
	newfile.close()
	
def rename_fastq_seq(fastq):
	dirc = os.path.dirname(fastq)
	name = os.path.basename(fastq)
	newfile = open(dirc+'/'+name.replace('fastq','renamed.fastq'), 'w')
	a=b=c=d=e=f=g=h=i = 0
	for line in open(fastq, 'r'):
		a += 1
		if line.startswith("@") and a%4 == 1:
			barcode = line.strip().split("=")[-1]
			newfile.write(line.replace('@','@'+barcode+'_'))
			b += 1
		else:
			newfile.write(line)
			c += 1	  
	newfile.close()
	print([a,b,c])

#rename_fasta_seq('/gs/project/wst-164-ab/anthony/Nanopore_data/olive_fly/TULIP_assembly/olivefly_illumina_seeds_100pb.fasta2')
#rename_fasta_seq2('/media/abayega/Padlock_DT/analysis/illumina_reads/olive_fly/bowtie/illumina_reads_mapping_supernova-pbjelly_once.fasta')
#rename_fasta_seq3()
#rename_fasta_seq4('/media/abayega/Padlock_DT/analysis/illumina_reads/olive_fly/bowtie/edited_illumina_reads_mapping_supernova-pbjelly_once.fasta')
#print(sys.argv[1])
#loop('/media/abayega/Padlock_DT/analysis/illumina_reads/olive_fly/bowtie/fasta_split')

rename_fastq_seq('/home/banthony/scratch/reads/nanopore/covid/B004/all_reads.fastq')