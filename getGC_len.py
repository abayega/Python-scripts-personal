#!/usr/bin/env python

import os, re

#This script takes a list of directories which contain many fasta files and loops over all files then returns sizes of DNA in bases and the GC content of the DNA
# export PATH=/gs/project/wst-164-ab/anthony/software/bin:$PATH
directories = ['/gs/project/wst-164-ab/anthony/test', '/gs/project/wst-164-ab/anthony/test2']
#dirc = ['/gs/project/wst-164-ab/anthony/test', '/gs/project/wst-164-ab/anthony/test2']

global sizes, GC_cont, readname, filename
#sizes = []
#GC_cont = []
#readname = []
#filename = []

def pursedir(directories):
    global sizes, GC_cont, readname, filename
    for i in directories:	
	getGC_len(i)

def getGC_len(directory):
	global sizes, GC_cont, readname, filename
	for x in directory:
	    sizes = []
	    GC_cont = []
	    readname = []
	    filename = []
	    print('You have declared ' + str(len(directory)) + ' directories to search, great!!!')
	    files = os.listdir(x)
	    print('You have ' + str(len(files)) + ' files in ' + str(os.path.basename(x)) + ', great!!!')
	    for i in range(len(files)):
		if files[i].endswith('.fasta'): # and 'r1' in files[i]:
			filename = filename + [files[i]]
			read = open(os.path.join(x, files[i]))
			fasta = read.readlines()
			for y in fasta:
				if re.match(r'^[G|T|C|A]', y):
					sizes = sizes + [len(y)]
					G=float(y.upper().count("G"))
					C=float(y.upper().count("C"))
					GC = (G+C)*100/len(y)					
					GC_cont = GC_cont + ["{0:.0f}".format(round(GC,0))]
					readname = readname + [fasta[fasta.index(y)-1].lstrip('>')]
					read.close()

	    data = open('{0}/%s_GCnSize.txt'.format(os.path.dirname(x)) % os.path.basename(x), 'a')
	    print('You have ' + str(len(sizes)) + ' files I will now write to disk')
	    for size in range(len(sizes)):
	    	data.write(str(filename[size]) + '\t' + str(readname[size]) + '\t' + str(sizes[size]) + '\t' + str(GC_cont[size]))
	    	data.write('\n')

	data.close()
	


#pursedir(directories)
getGC_len(directories)
#print(sizes)
	
