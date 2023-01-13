#!/usr/bin/env python
import re, sys, os, progressbar

# vrlpr.py: find overlapping genes in a GTF file
# Usage: python vrlpr.py < genes.gtf > overlaps.txt

def overlap(range1, range2):
	return range1[0] == range2[0] and range1[2] >= range2[1] and range1[1] <= range2[2]

def overlap_genes(gtf):
	dirc = os.path.dirname(gtf)
	name = os.path.basename(gtf)
	nfile1 = open(dirc+'/'+name+'_overlapping_genes', 'w')
	#nfile2 = open(dirc+'/'+'gtf_new2', 'w')
	a=b=c=d=e=f=g= 0
	genes = {}
	length = int(os.popen('wc -l %s' %(gtf)).read().split()[0])
	print('\nWe are looping over ' + str(length) + ' lines\n')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(gtf, 'r'):
		#print(line)
		if line == "" or line.startswith("#") or line.split('\t')[2] != 'gene':
			continue
		line = line.rstrip()
		fields = line.split("\t")
		seqid = fields[0]
		start = fields[3] #int(fields[3])
		end = fields[4] #int(fields[4])
		attributes = fields[8]

		matches = re.match('gene_id "([^"]+)"', attributes)
		if matches:
			geneid = matches.group(1)
			if geneid not in genes:
				genes[geneid] = [seqid, start, end, geneid]
			else:
				assert seqid == genes[geneid][0]
				if start < genes[geneid][1]:
					genes[geneid][1] = start
				if end > genes[geneid][2]:
					genes[geneid][2] = end
		'''
		for i in range(len(genes)):
			for j in range(i+1, len(genes)):
				geneid_i = genes.keys()[i]
				geneid_j = genes.keys()[j]
				range_i = genes[geneid_i]
				range_j = genes[geneid_j]
				if overlap(range_i, range_j):
					nfile1.write("%s\t%s" %(geneid_i, geneid_j))
		'''			
		d+=1
		bar.update(d)
	bar.finish()

	d=0
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(gtf, 'r'):
		if line == "" or line.startswith("#") or line.split('\t')[2] != 'gene':
			continue
		line = line.rstrip()
		fields = line.split("\t")
		seqid = fields[0]
		start = fields[3]
		end = fields[4]
		attributes = fields[8]

		matches = re.match('gene_id "([^"]+)"', attributes)
		if matches:
			geneid = matches.group(1)
			for x,y in genes.items():
				start2 = y[1]
				end2 = y[2]
				if int(start) >= int(start2) and int(start) <= int(end2) or int(end) >= int(start2) and int(end) <= int(end2):
					if geneid != x and seqid == y[0]:
						nfile1.write('\t'.join([seqid, start, end, geneid, 'n','\t'.join(y)])+'\n')
		d+=1
		bar.update(d)
	bar.finish()
	nfile1.close()

def overlap_introns_exons(intron, exons):
	dirc = os.path.dirname(intron)
	name = os.path.basename(intron)
	nfile1 = open(dirc+'/'+'introns_exons_overlapping_genes', 'w')
	#nfile2 = open(dirc+'/'+'gtf_new2', 'w')
	a=b=c=d=e=f=g= 0
	genes = {}
	length = int(os.popen('wc -l %s' %(intron)).read().split()[0])
	print('\nWe are looping over ' + str(length) + ' lines\n')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(intron, 'r'):
		#print(line)
		if line == "" or line.startswith("#"): # or line.split('\t')[2] != 'gene':
			continue
		line = line.rstrip()
		fields = line.split("\t")
		seqid = fields[0]
		start = fields[1] #int(fields[3])
		end = fields[2] #int(fields[4])
		#attributes = fields[8]

		#matches = re.match('gene_id "([^"]+)"', attributes)
		matches = fields[3].startswith('gene')
		if matches:
			geneid = fields[3] #matches.group(1)
			if geneid not in genes:
				genes[geneid] = [seqid, start, end, geneid]
			else:
				assert seqid == genes[geneid][0]
				if start < genes[geneid][1]:
					genes[geneid][1] = start
				if end > genes[geneid][2]:
					genes[geneid][2] = end		
		d+=1
		bar.update(d)
	bar.finish()
	length = int(os.popen('wc -l %s' %(exons)).read().split()[0])
	d=0
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(exons, 'r'):
		if line == "" or line.startswith("#"): # or line.split('\t')[2] != 'gene':
			continue
		line = line.rstrip()
		fields = line.split("\t")
		seqid = fields[0]
		start = fields[1]
		end = fields[2]

		matches = fields[3].startswith('gene')
		if matches:
			geneid = fields[3]
			for x,y in genes.items():
				start2 = y[1]
				end2 = y[2]
				if int(start) >= int(start2) and int(start) <= int(end2) or int(end) >= int(start2) and int(end) <= int(end2):
					if geneid != x and seqid == y[0]:
						nfile1.write('\t'.join([seqid, start, end, geneid, 'n','\t'.join(y)])+'\n')
		d+=1
		bar.update(d)
	bar.finish()
	nfile1.close()
					
					
dr = '/home/banthony/scratch/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/rna_velocity/hamed/annotation'
#overlap_genes(dr + '/' + 'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds.gff3.mrna.gencode.gtf.2')
#overlap_genes(dr + '/' + 'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds.gff3.mrna.gencode.gtf.collapsed.2')
overlap_introns_exons(dr + '/' + 'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds.gff3.mrna.gencode.gtf.collapsed.genePred.sorted.introns.bed12',
					  dr + '/' + 'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds.gff3.mrna.gencode.gtf.collapsed.genePred.sorted.bed12')
				
