#!/usr/bin/env python3

import re, os, sys, shelve #, patternCount5
import progressbar

def fasta_to_dict(fasta):
	#This function takes a fastq file with single line sequence and returns a dictionary with sequence name and corresponding read
	fastaNameSeq = {}	
	openfasta = open(fasta, 'r')
	
	x = 0
	
	length = int(os.popen('wc -l %s' %(fasta)).read().split()[0])
	
	file1=open(fasta,'r')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length:
		a=file1.readline()
		b=file1.readline()
		
		fastaNameSeq.setdefault(a, b)		
		x+=2
		bar.update(x)
	bar.finish()
	
	#print(fastqNameSeq)
	return(fastaNameSeq)

def fastq_to_dict(fastq):
	#This function takes a fastq file with single line sequence and returns a dictionary with sequence name and corresponding read
	fastqNameSeq = {}	
	openfastq = open(fastq, 'r')
	
	x = 0
	
	length = int(os.popen('wc -l %s' %(fastq)).read().split()[0])
	file1=open(fastq,'r')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length:
		a=file1.readline()
		b=file1.readline()
		c=file1.readline()
		d=file1.readline()
		
		fastqNameSeq.setdefault(a, b+c+d)		
		x+=4
		bar.update(x)
	bar.finish()
	
	#print(fastqNameSeq)
	return(fastqNameSeq)

def look_up_given_reads_fasta(file, fasta): #Use this for speed
	#This script takes a fasta file and a file with the names of reas u want to take from the fasta. It will create a dictionary for the fasta seqs and search read names against the fasta name. 
	#If a name in the file is found in the fasta as well, the corresponding sequence and the name are returned. Pay attention to the \n and > and fasta format, it has to be single line for seqs
	lopenfasta = fasta_to_dict(fasta)
	openfile=open(file,'r')
	lopenfile = openfile.readlines()
	l = len(lopenfile)
	path = os.path.dirname(file)
	name = os.path.basename(file)
	outq=open(path+'/'+name+'_2.fasta','w')	
	
	n = 0
	x = 0
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for read in lopenfile:		
		if read in lopenfasta.keys():
			outq.write(read+lopenfasta[read])
			x += 1
		n+=1
		bar.update(n)
	bar.finish()
	outq.close()
	print('Discovered ' + str(x) + ' sequences' + ' which is %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*x/l),1))))	

def look_up_given_reads_fastq(file, fastq): #Use this for speed
	#This script takes a fasta file and a file with the names of reas u want to take from the fasta. It will create a dictionary for the fasta seqs and search read names against the fasta name. 
	#If a name in the file is found in the fasta as well, the corresponding sequence and the name are returned. Pay attention to the \n and > and fasta format, it has to be single line for seqs
	lopenfastq = fastq_to_dict(fastq)
	openfile=open(file,'r')
	lopenfile = openfile.readlines()
	l = len(lopenfile)
	path = os.path.dirname(file)
	name = os.path.basename(file)
	outq=open(path+'/'+name+'_2.fastq','w')	
	
	n = 0
	x = 0
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for read in lopenfile:		
		if read.replace('>','@') in lopenfastq.keys():
			outq.write(read.replace('>','@')+lopenfastq[read.replace('>','@')])
			x += 1
		n+=1
		bar.update(n)
	bar.finish()
	outq.close()
	print('Discovered ' + str(x) + ' sequences' + ' which is %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*x/l),1))))
	
def look_up_given_reads_lengths(file, fastx, file_format):
	openfile=open(file,'r')
	lopenfile = openfile.readlines()
	l = len(lopenfile)
	path = os.path.dirname(file)
	name = os.path.basename(file)
	outx=open(path+'/'+name+'_2.txt','w')
	
	n = 0
	x = 0
	
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	if file_format == 'fasta':
		lopenfasta = fasta_to_dict(fastx)
		fastxq = 0
	elif file_format == 'fastq':
		lopenfastq = fastq_to_dict(fastx)
		fastxq = 1
		
	for read in lopenfile:		
		if fastxq == 1 and read.replace('>','@') in lopenfastq.keys():
			outx.write(read.replace('>','@').strip()+'\t'+str(len(lopenfastq[read.replace('>','@')].split('\n')[0]))+'\n')
			x += 1
			
		elif fastxq == 0 and read in lopenfasta.keys():
			outx.write(read.strip()+'\t'+str(len(lopenfastq[read])))
			x += 1
		n+=1
		bar.update(n)		
	bar.finish()
	outx.close()
	print('Discovered ' + str(x) + ' sequences' + ' which is %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*x/l),1))))
	
def return_read_names_lengths(fastx, file_format):
	path = os.path.dirname(fastx)
	name = os.path.basename(fastx).replace('.'+file_format, '')
	outx=open(path+'/'+name+'_readlengths.txt','w')
	
	n = 0
	xi = 0
	
	if file_format == 'fasta':
		lopenfasta = fasta_to_dict(fastx)
		l = len(lopenfasta.keys())
		fastxq = 0
		bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		bar.start()
		for x,y in lopenfasta.items():
			outx.write(x.strip()+'\t'+str(len(y.strip()))+'\n')
			xi += 1
			n+=1
			bar.update(n)
	elif file_format == 'fastq':
		lopenfastq = fastq_to_dict(fastx)
		l = len(lopenfasta.keys())
		fastxq = 1
		bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		bar.start()
		for x,y in lopenfasta.items():
			outx.write(x.strip()+'\t'+str(len(y.split('\n')[0].strip()))+'\n')
			xi += 1
			n+=1
			bar.update(n)		
	bar.finish()
	outx.close()
	print('Discovered ' + str(xi) + ' sequences' + ' which is %s percent of all sequences in the file' %(str("{0:.1f}".format(round(100*xi/l),1))))

	
dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/minimap2_all_analyse'
#look_up_given_reads_fasta('/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/ercc/gmap/mapped_reads',\
#					'/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/ercc/gmap/Bo.E.5H_all_rnaseq_combined_pass.trimmed_stranded2.fasta')
#look_up_given_reads_fastq(dr+'/'+'mapped_reads',dr+'/'+'Bo.E.5H_all_rnaseq_combined_pass.trimmed_stranded2_5522273-7363016.fastq')
#look_up_given_reads_lengths('sequences.mpa-format.labels.sorted4','/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/illumina/olivefly/kraken_male_genome/14.txt')
return_read_names_lengths(dr+'/'+'best.sorted2_sorted_collapsed_genes.fasta', 'fasta')