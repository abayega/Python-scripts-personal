#!/usr/bin/env python3

import re, os, progressbar, sys

def convertFasta(fasta):
	#newfile_name = os.path.join(os.path.dirname(fasta), '_sl.'.join([os.path.basename(fasta).replace('.'+fileformat,''),fileformat]))
	base = os.path.basename(fasta)
	if base.endswith('fasta'):
		newfile_name = os.path.dirname(fasta) + '/' + os.path.basename(fasta).replace('.fasta','_sl.fasta')
		newfile = open(newfile_name, 'w')
	elif base.endswith('fa'):
		newfile_name = os.path.dirname(fasta) + '/' + os.path.basename(fasta).replace('.fa','_sl.fa')
		newfile = open(newfile_name, 'w')
	else:
		#print('Fasta file must either end with fasta or fa, script has crashed')
		newfile_name = os.path.dirname(fasta) + '/' + os.path.basename(fasta)+'_sl'
		newfile = open(newfile_name, 'w')
	#if not os.path.exists(newfile_name):
		#newfile = open(newfile_name, 'w')

	reads = int(os.popen('grep ">" %s | wc -l' %(fasta)).read().split()[0])
	print('\nDetected ' + str(reads) + ' reads\n')
	length = int(os.popen('wc -l %s' %(fasta)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines\n')
	file1=open(fasta,'r')
	k=0
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while k < length:
		a=file1.readline()
		if not a.startswith('>'):
			newfile.write(a.strip())
		elif k==0 and a.startswith('>'):
			newfile.write(a)
		else:
			newfile.write('\n'+a)
		k+=1
		bar.update(k)
	bar.finish()
	newfile.close()
	return(newfile_name)
	
def check_fasta_fmt(infile):
	path = os.path.dirname(infile)
	#Check first
	length1 = int(os.popen('head -40000 %s | grep ">" | wc -l' %(infile)).read().split()[0])
	length2 = int(os.popen('head -40000 %s | grep -v ">" | wc -l' %(infile)).read().split()[0])
	if length1 == length2:
		print('\nFasta format is good...\n')
		newfile_name = infile
	else:
		print('\nFasta format is bad, creating a single line format, delete file if not needed...')
		newfile_name=convertFasta(infile)
		#infile = path + '/' + os.path.basename(infile).replace('.fasta','_sl.fasta')
		#infile = newfile_name
	return(str(newfile_name))
