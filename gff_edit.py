#!/usr/bin/env python3

import os, re
import progressbar
			
def gff3_to_gtf(gff3):
	newfile = open(os.path.join(os.path.dirname(gff3), '_'.join([os.path.basename(gff3).replace('.gff3',''),'edited.gff3'])), 'w')
	gff3_file = open(gff3, 'r')
	lgff3_file = gff3_file.readlines()
	for i in lgff3_file:
		if 'gene' in i or 'mRNA' in i:
			gene2 = []
			split = i.split('\t')
			gene2 = split[:8]
			ID = split[8].split(';')
			j = []
			for m in ID:
				if m.startswith('ID='):
					j += [m]
				elif m.startswith('Name='):
					j += [m]
				elif m.startswith('Parent='):
					j += [m]
				
			k = ';'.join(j)
			gene2+=[k]
			newfile.write('\t'.join(gene2))
			newfile.write('\n')
		else:
			newfile.write(i)
			

gff3_to_gtf('/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/capitata/rna_seq/CCAP.Models.gff3')
