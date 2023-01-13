#!/usr/bin/env python3

import os, re
import progressbar

def collapse_genes(bed,fileformat):
	print('Starting to collaspse genes')
	obed = open(bed)
	rbed = obed.readlines()	
	newfile = open(os.path.join(os.path.dirname(bed), '_2.'.join([os.path.basename(bed).replace('.'+fileformat,''),fileformat])),'w')
	k = 0
	genes = []
	pre_gene = []
	features = []
	length = []
	lrbed = len(rbed)
	bar = progressbar.ProgressBar(maxval=len(rbed), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in rbed:
		#print(line.split('\t'))
		index = rbed.index(line)
		start = line.split('\t')[1]
		end = line.split('\t')[2]
		if rbed.index(line) < lrbed-1:
			lnext = rbed[index+1]
			if line.split('\t')[0] == lnext.split('\t')[0] and line.split('\t')[5] == lnext.split('\t')[5]:			
				lstart = lnext.split('\t')[1]
				lend = lnext.split('\t')[2]
				if int(lstart) in range(int(start),int(end)+1):
					pre_gene += [line]
				else:
					pre_gene += [line]
					for seq in pre_gene:
						start1 = seq.split('\t')[1]
						end1 = seq.split('\t')[2]
						length += [int(end1) - int(start1)]
					longest = pre_gene[length.index(max(length))]
					genes += [longest]
					pre_gene = []
					length = []
			elif rbed.index(line) != 0 and line.split('\t')[0] != lnext.split('\t')[0] and line.split('\t')[0] == rbed[index-1].split('\t')[0]:
				lprev = rbed[index-1]
				lpstart = lprev.split('\t')[1]
				lpend = lprev.split('\t')[2]
				if int(line.split('\t')[1]) in range(int(lpstart),int(lpend)+1):
					pre_gene += [line]
					for seq in pre_gene:
						start1 = seq.split('\t')[1]
						end1 = seq.split('\t')[2]
						length += [int(end1) - int(start1)]
					longest = pre_gene[length.index(max(length))]
					genes += [longest]
					pre_gene = []
					length = []
				else:
					genes += [line]
			elif rbed.index(line) == 0 and line.split('\t')[0] != lnext.split('\t')[0]:
				genes += [line]
			elif rbed.index(line) != 0 and line.split('\t')[0] != lnext.split('\t')[0] and line.split('\t')[0] != rbed[index-1].split('\t')[0]:
				genes += [line]
			#elif rbed.index(line) != 0 and line.split('\t')[0] != lnext.split('\t')[0] and line.split('\t')[0] == rbed[index-1].split('\t')[0]:
			#	genes += [line]
				
		elif line.split('\t')[0] == rbed[index - 1].split('\t')[0]:
			pre_gene += [line]
			for seq in pre_gene:
				start1 = seq.split('\t')[1]
				end1 = seq.split('\t')[2]
				length += [int(end1) - int(start1)]
			longest = pre_gene[length.index(max(length))]
			genes += [longest]
			pre_gene = []
			length = []
		else:
			v = 0
			genes += [line]
		k+=1
		bar.update(k)
	bar.finish()
	print('Now writing results...')
	for i in genes:
		newfile.write(i)
	newfile.close
	
dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/minimap2_all_analyse'

collapse_genes(dr+'/'+'best.sorted2_sorted_1-308219.bed','bed')