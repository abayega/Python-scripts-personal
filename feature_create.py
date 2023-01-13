#!/usr/bin/env python3

import re, os, progressbar

def feature_create(fasta):
	#Script takes a fasta file then creates a bed that has the read name, 0, and length of sequence all tab delimited
	ofasta = open(fasta, 'r')
	lenfasta = int(os.popen('wc -l %s' %(fasta)).read().split()[0])
	out = open(fasta+'_bed', 'w')
	x = 0
	
	bar = progressbar.ProgressBar(maxval=lenfasta, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < lenfasta:
		a=ofasta.readline().strip()
		b=ofasta.readline().strip()
		
		out.write(a.strip('>').split(' ')[0] + '\t' + '0' + '\t' + str(len(b)) + '\n')
		
		x += 2
		
		bar.update(x)
	bar.finish()
	out.close()
	
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/combined/rnaseq_saturation'

feature_create(dr + '/' + 'ERCC92_ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_gffread_sorted_edited_over37_transcripts_sl2.fasta')
	
		
		
	
	