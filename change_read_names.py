#!/usr/bin/env python3

import re, os, progressbar

def change_read_names(reads_fasta,gpd):
	dirc = os.path.dirname(reads_fasta)
	nfile1 = open(reads_fasta+'_2', 'w')
	nfile2 = open(gpd+'_2', 'w')
	
	
	x=c=d=e=f=g=h=i = 0
	
	length = int(os.popen('wc -l %s' %(reads_fasta)).read().split()[0])
	print('We are looping over ' + str(length/2) + ' reads')
	file1=open(reads_fasta,'r')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length:
		a=file1.readline()
		b=file1.readline()
		nfile1.write('_'.join(a.strip().split('.'))+'\n'+b)
		c += 1			
		x+=2
		bar.update(x)
	bar.finish()
	if length/2 != sum([c,d,e,f,g,h,i]):
		print('The total lines in file does not equate number of reads distributed\n' + str(length/2) + ' ' + str(sum([c,d,e,f,g,h,i])))
	else:
		print('The total lines in file is ' + str(length/2))
	
	file2=open(gpd,'r')
	length = int(os.popen('wc -l %s' %(gpd)).read().split()[0])
	print('We are looping over ' + str(length) + ' gpd lines')
	x = 0
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length:
		a=file2.readline()
		a = a.strip().split()
		nfile2.write('\t'.join([a[0],'_'.join(a[1].split('.'))])+ '\t'.join(a[2:])+b)
		d += 1			
		x+=1
		bar.update(x)
	bar.finish()
	if length != d:
		print('The total lines in file does not equate number of reads distributed\n' + str(length/2) + ' ' + str(sum([c,d,e,f,g,h,i])))
	else:
		print('The total lines in file is ' + str(length))