#!/usr/bin/env python

import os, re

#This script takes a list of directories which contain many fasta files and loops over all files then returns sizes of DNA in bases

directories = ['/gs/project/wst-164-ab/anthony/Nanopore_data/C006_01_5/nanook/fasta/pass/Template']

sizes = []

for x in directories:
    print('You have declared ' + str(len(directories)) + ' directories to search, great!!!')
    files = os.listdir(x)
    print('You have ' + str(len(files)) + ' files in ' + str(os.path.basename(x)) + ', great!!!')
    for i in range(len(files)):
	if files[i].endswith('.fasta') and 'r1' in files[i]:
		read = open(os.path.join(x, files[i]))
		fasta = read.readlines()
		for y in fasta:
			if re.match(r'^[G|T|C|A]', y):
				sizes = sizes + [len(y)]
				read.close()

data = open('r1Template_fastasize.txt', 'a')
for size in sizes:
    data.write(str(size))
    data.write('\n')

data.close()
print(len(sizes))
#print(sizes)
	
