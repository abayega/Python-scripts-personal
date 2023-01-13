#!/usr/bin/env python

#This code takes a file containing mulitple lines with each line containing a number. It bins the numbers and prints to std.out the sum of the bins

import os

r1k = 0
r1-5k = 0
r5-10k = 0
r10-20k = 0
r20-30k = 0
r30-50k = 0
r>50k = 0


read = open(doc)
fasta = read.readlines()
for x in fasta:
    if int(x) <= 1000:
	r1k = r1k + int(x)
    elif int(x) > 1000 and int(x) <= 5000:
	r1-5k = r1-5k + int(x):
    elif int(x) > 5000 and int(x) <= 10000:
	r5-10k = r5-10k + int(x):
    elif int(x) > 10000 and int(x) <= 20000:
	r10-20k = r10-20k + int(x):
    elif int(x) > 20000 and int(x) <= 30000:
	r20-30k = r20-30k + int(x):
    elif int(x) > 30000 and int(x) <= 50000:
	r30-50k = r30-50k + int(x):
    elif int(x) > 50000:
	r>50k = r>50k + int(x):

print(r1h, r1-5k, r5-10k, r10-20k, r20-30k, r30-50k, r>50k)
	
