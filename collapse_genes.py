#!/usr/bin/env python3

import os, re
import progressbar
	
def collapse_genes(bed,fileformat):
	global pre_gene, findftr, ofindftr, features, isoforms
	print('\nStarting to collaspse genes\n')
	obed = open(bed)
	rbed = obed.readlines()	
	newfile = open(os.path.join(os.path.dirname(bed), '_collapsed_genes.'.join([os.path.basename(bed).replace('.'+fileformat,''),fileformat])),'w')
	findftr = os.path.join(os.path.dirname(bed), '_delete.'.join([os.path.basename(bed).replace('.'+fileformat,''),fileformat]))
	error_lines = open(os.path.join(os.path.dirname(bed), '_errors.'.join([os.path.basename(bed).replace('.'+fileformat,''),fileformat])),'w')
	#ofindftr = open(findftr,'w')
	k = 0
	genes = []
	#isoforms = []
	pre_gene = []
	length = []
	
	lrbed = len(rbed)
	bar = progressbar.ProgressBar(maxval=len(rbed), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in rbed:
		try:
			features = []
			isoforms = []
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
						a = pre_gene[0].split('\t')[1]
						b = 0
						for seq in pre_gene:
							start1 = seq.split('\t')[1]
							end1 = seq.split('\t')[2]
							length += [int(end1) - int(start1)]
							b = max(int(end1),b)
						longest = pre_gene[length.index(max(length))]
						gene_features = get_gene_features(bed)
						got_isoforms = get_isoforms(bed)
						itms = longest.split('\t')
						itms2 = gene_features.strip().split('\t')
						genes += [itms[0]+'\t'+a+'\t'+str(b)+'\t'+itms[3]+'\t'+itms[4]+'\t'+itms[5]+'\t'+a+'\t'+str(b)+'\t'+itms[8]+'\t'+itms2[2]+'\t'+itms2[0]+'\t'+itms2[1]+'\t'+got_isoforms+'\n']
						pre_gene = []
						length = []
				elif rbed.index(line) != 0 and line.split('\t')[0] != lnext.split('\t')[0] and line.split('\t')[0] == rbed[index-1].split('\t')[0]:
					lprev = rbed[index-1]
					lpstart = lprev.split('\t')[1]
					lpend = lprev.split('\t')[2]
					if int(line.split('\t')[1]) in range(int(lpstart),int(lpend)+1):
						pre_gene += [line]
						a = pre_gene[0].split('\t')[1]
						b = 0
						for seq in pre_gene:
							start1 = seq.split('\t')[1]
							end1 = seq.split('\t')[2]
							length += [int(end1) - int(start1)]
							b = max(int(end1),b)
						longest = pre_gene[length.index(max(length))]
						gene_features = get_gene_features(bed)
						got_isoforms = get_isoforms(bed)
						itms = longest.split('\t')
						itms2 = gene_features.strip().split('\t')
						#print(gene_features)
						#print(itms2[2]) #+'\t'+itms2[0]+'\t'+itms2[1])
						genes += [itms[0]+'\t'+a+'\t'+str(b)+'\t'+itms[3]+'\t'+itms[4]+'\t'+itms[5]+'\t'+a+'\t'+str(b)+'\t'+itms[8]+'\t'+itms2[2]+'\t'+itms2[0]+'\t'+itms2[1]+'\t'+got_isoforms+'\n']
						pre_gene = []
						length = []
					else:
						a = line.split('\t')[1]
						b = line.split('\t')[2]
						pre_gene += [line]
						gene_features = get_gene_features(bed)
						got_isoforms = get_isoforms(bed)
						itms = line.split('\t')
						itms2 = gene_features.strip().split('\t')
						genes += [itms[0]+'\t'+a+'\t'+str(b)+'\t'+itms[3]+'\t'+itms[4]+'\t'+itms[5]+'\t'+a+'\t'+str(b)+'\t'+itms[8]+'\t'+itms2[2]+'\t'+itms2[0]+'\t'+itms2[1]+'\t'+got_isoforms+'\n']
						pre_gene = []
						length = []
				elif rbed.index(line) == 0 and line.split('\t')[0] != lnext.split('\t')[0]:
					a = line.split('\t')[1]
					b = line.split('\t')[2]
					pre_gene += [line]
					gene_features = get_gene_features(bed)
					got_isoforms = get_isoforms(bed)
					itms = line.split('\t')
					itms2 = gene_features.strip().split('\t')
					genes += [itms[0]+'\t'+a+'\t'+str(b)+'\t'+itms[3]+'\t'+itms[4]+'\t'+itms[5]+'\t'+a+'\t'+str(b)+'\t'+itms[8]+'\t'+itms2[2]+'\t'+itms2[0]+'\t'+itms2[1]+'\t'+got_isoforms+'\n']
					pre_gene = []
					length = []
				elif rbed.index(line) != 0 and line.split('\t')[0] != lnext.split('\t')[0] and line.split('\t')[0] != rbed[index-1].split('\t')[0]:
					a = line.split('\t')[1]
					b = line.split('\t')[2]
					pre_gene += [line]
					gene_features = get_gene_features(bed)
					got_isoforms = get_isoforms(bed)
					itms = line.split('\t')
					itms2 = gene_features.strip().split('\t')
					genes += [itms[0]+'\t'+a+'\t'+str(b)+'\t'+itms[3]+'\t'+itms[4]+'\t'+itms[5]+'\t'+a+'\t'+str(b)+'\t'+itms[8]+'\t'+itms2[2]+'\t'+itms2[0]+'\t'+itms2[1]+'\t'+got_isoforms+'\n']
					pre_gene = []
					length = []
					b=0
			elif line.split('\t')[0] == rbed[index - 1].split('\t')[0]:
				#print(1)
				pre_gene += [line]
				a = pre_gene[0].split('\t')[1]
				b = 0
				for seq in pre_gene:
					start1 = seq.split('\t')[1]
					end1 = seq.split('\t')[2]
					length += [int(end1) - int(start1)]
					b = max(int(end1),b)
				longest = pre_gene[length.index(max(length))]
				gene_features = get_gene_features(bed)
				got_isoforms = get_isoforms(bed)
				itms = longest.split('\t')
				itms2 = gene_features.strip().split('\t')
				genes += [itms[0]+'\t'+a+'\t'+str(b)+'\t'+itms[3]+'\t'+itms[4]+'\t'+itms[5]+'\t'+a+'\t'+str(b)+'\t'+itms[8]+'\t'+itms2[2]+'\t'+itms2[0]+'\t'+itms2[1]+'\t'+got_isoforms+'\n']
				pre_gene = []
				length = []
			else:

				a = line.split('\t')[1]
				b = line.split('\t')[2]
				pre_gene += [line]
				gene_features = get_gene_features(bed)
				got_isoforms = get_isoforms(bed)
				itms = longest.split('\t')
				itms2 = gene_features.strip().split('\t')
				genes += [itms[0]+'\t'+a+'\t'+str(b)+'\t'+itms[3]+'\t'+itms[4]+'\t'+itms[5]+'\t'+a+'\t'+str(b)+'\t'+itms[8]+'\t'+itms2[2]+'\t'+itms2[0]+'\t'+itms2[1]+'\t'+got_isoforms+'\n']
				pre_gene = []
				length = []
		except:
			error_lines.write(line.strip() + '\n')
		k+=1
		bar.update(k)
	bar.finish()
	print('\nNow writing results...\n')
	for i in genes:
		newfile.write(i)
	newfile.close 
	
def get_gene_features(bed):
	global pre_gene, findftr, ofindftr, features, isoforms
	ofindftr = open(findftr,'w')
	for line in pre_gene:
		noftr = line.split('\t')[9]
		lenftr = line.split('\t')[10].strip(',').split(',')
		startftr = line.split('\t')[11].strip(',').split(',')
		for x in range(len(lenftr)):
			#ofindftr.write(str(startftr[x]) + '\t' + str(str(startftr[x])+lenftr[x]))
			ofindftr.write(startftr[x] + '\t' + str(int(startftr[x])+int(lenftr[x]))+'\n')
	ofindftr.close()
	os.system('sort -k1,1n %s > %s_sorted' %(findftr,findftr))
	sortedftr = open(findftr+'_sorted', 'r')
	rsortedftr = sortedftr.readlines()
	pre_features = []
	feature = ':'.join(rsortedftr[0].split('\t')).strip()
	for eline in rsortedftr:
		start2 = eline.split('\t')[0]
		end2 = eline.split('\t')[1].strip()
		j = int(feature.split(':')[0])
		k = int(feature.split(':')[1])
		#print([k])
		if int(start2) in range(j,k+1) and int(end2) in range(j,k+1):
			uuu = 0 # forget about this feature
		elif int(start2) in range(j,k+1) and not int(end2) in range(j,k+1):
			feature = str(j)+':'+ end2
			feature = feature.strip()
		elif not int(start2) in range(j,k+1) and not int(end2) in range(j,k+1):
			features += [feature]
			feature = ':'.join(eline.split('\t'))
			feature = feature.strip()
			
		if rsortedftr.index(eline) == len(rsortedftr) -1:
			features += [feature]
			
	no_ftr = len(features)
	size_ftr = []
	start_ftr = []
	for ftr in features:
		start_ftr += [ftr.split(':')[0]]
		size_ftr += [str(int(ftr.split(':')[1])-int(ftr.split(':')[0]))]
	start_ftr += ['']
	start_ftr = ','.join(start_ftr)
	size_ftr += ['']
	size_ftr = ','.join(size_ftr)
	
	return (size_ftr+'\t'+start_ftr+'\t'+str(no_ftr))	
	#return ','.join(features)
	
def get_isoforms(bed):
	global pre_gene, findftr, ofindftr, features, isoforms
	ofindftr = open(findftr,'w')
	for line in pre_gene:
		noftr = line.split('\t')[9]
		lenftr = line.split('\t')[10].strip(',').split(',')
		startftr = line.split('\t')[11].strip(',').split(',')
		ofindftr.write(noftr + '\t')
		for x in range(len(lenftr)):
			#ofindftr.write(str(startftr[x]) + '\t' + str(str(startftr[x])+lenftr[x]))
			ofindftr.write(startftr[x] + ':' + str(int(startftr[x])+int(lenftr[x]))+'\t')
		ofindftr.write('\n')		
	ofindftr.close()
	os.system('sort -u -k1,1n -k2,2n %s > %s_sorted_iso' %(findftr,findftr))
	sortedftr = open(findftr+'_sorted_iso', 'r')
	rsortedftr = sortedftr.readlines()
	feature = rsortedftr[0].strip().split('\t')
	for eline in rsortedftr:
		eline1 = eline.strip().split('\t')
		if feature[0] == eline1[0]:			
			for i in range(int(feature[0])):
				if int(eline1[i+1].split(':')[0]) in range(int(feature[i+1].split(':')[0]) - 5, int(feature[i+1].split(':')[0]) + 5) and int(eline1[i+1].split(':')[1]) in range(int(feature[i+1].split(':')[1]) - 5, int(feature[i+1].split(':')[1]) + 5):
					moveon = 1 # we move on to next exon
				#elif not int(eline1[i+1].split(':')[0]) in range(int(feature[i+1].split(':')[0]) - 5, int(feature[i+1].split(':')[0]) + 5):
				else:
					size_ftr = []
					start_ftr = []
					for ftr in feature[1:]:
						start_ftr += [ftr.split(':')[0]]
						size_ftr += [str(int(ftr.split(':')[1])-int(ftr.split(':')[0]))]
					start_ftr += ['']
					start_ftr = ','.join(start_ftr)
					size_ftr += ['']
					size_ftr = ','.join(size_ftr)
					isoforms += [feature[0]+'-'+size_ftr+':'+start_ftr+'\t']
					#print(feature[0]+'-'+size_ftr+'-'+start_ftr+'\t')
					feature = eline1
					
		else:
			size_ftr = []
			start_ftr = []
			for ftr in feature[1:]:
				start_ftr += [ftr.split(':')[0]]
				size_ftr += [str(int(ftr.split(':')[1])-int(ftr.split(':')[0]))]
			start_ftr += ['']
			start_ftr = ','.join(start_ftr)
			size_ftr += ['']
			size_ftr = ','.join(size_ftr)
			isoforms += [feature[0]+'-'+size_ftr+':'+start_ftr+'\t']
			#print(feature[0]+'-'+size_ftr+'-'+start_ftr+'\t')
			feature = eline1
			
		if rsortedftr.index(eline) == len(rsortedftr) -1:
			size_ftr = []
			start_ftr = []
			for ftr in feature[1:]:
				start_ftr += [ftr.split(':')[0]]
				size_ftr += [str(int(ftr.split(':')[1])-int(ftr.split(':')[0]))]
			start_ftr += ['']
			start_ftr = ','.join(start_ftr)
			size_ftr += ['']
			size_ftr = ','.join(size_ftr)
			isoforms += [feature[0]+'-'+size_ftr+':'+start_ftr+'\t']
			
			
	return ('\t'.join(isoforms))

def parse_bed(bed,fileformat):
	#This script parses the bed file and removes bad lines, basically cleaning it up to prepare for collapsing genes
	print('Starting to parse bed file')
	obed = open(bed)
	rbed = obed.readlines()
	rbed = rbed[1:]
	newbed = open(os.path.join(os.path.dirname(bed), '_edited.'.join([os.path.basename(bed).replace('.'+fileformat,''),fileformat])),'w')
	newbed2 = open(os.path.join(os.path.dirname(bed), '_edited_bad.'.join([os.path.basename(bed).replace('.'+fileformat,''),fileformat])),'w')
	k = 0
	m = 0
	n = 0
	o = 0
	lrbed = len(rbed)
	bar = progressbar.ProgressBar(maxval=lrbed, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in rbed:
		intron = []
		linex = line.split('\t')		
		if int(linex[9]) > 1:
			size = linex[10].split(',')
			start = linex[11].split(',')
			if int(size[0]) < 50 and int(start[1]) > 5000:
				newbed2.write(line)
				m += 1
			elif int(size[0]) <= 10:
				newbed2.write(line)
				m += 1
			elif int(linex[9]) >= 2 and int(size[len(size)-2]) < 50 and int(start[len(start)-2]) - (int(size[len(size)-3]) + int(start[len(start)-3])) > 5000:
				#print(int(start[len(start)-2]))
				newbed2.write(line)
				m += 1
			elif int(linex[9]) > 3 and int(size[len(size)-2]) <= 10:
				newbed2.write(line)
				m += 1
			elif int(start[len(start)-2]) - int(start[0]) > 50000:
				newbed2.write(line)
				m += 1
			else:
				newbed.write(line)
				n += 1
		else:
			newbed.write(line)
			o += 1
		k+=1
		bar.update(k)
	bar.finish()
	newbed.close()
	newbed2.close()
	fname = os.path.join(os.path.dirname(bed), '_edited.'.join([os.path.basename(bed).replace('.'+fileformat,''),fileformat]))
	newfname = os.path.join(os.path.dirname(fname), '_resorted.'.join([os.path.basename(fname).replace('.'+fileformat,''),fileformat]))
	os.system('sort -k6,6 -k1,1 -k2,2n %s > %s' %(fname,newfname))
	print('There were ' + str(m) + ' bad lines which is ' + str("{0:.2f}".format(round(m/lrbed * 100,2))) + '% of the total lines.' + '\n' \
		  + 'There were ' + str(n) + ' good lines which is ' + str("{0:.0f}".format(round(n/lrbed * 100,0))) + '% of the total lines')
	print('There were ' + str(o) + ' lines with 1 feature which is ' + str("{0:.0f}".format(round(o/lrbed * 100,0))) + '% of the total lines.')
	print('Total lines in file are ' + str(lrbed))
	
def gene_mapping(bed, fileformat):
	obed = open(bed)
	rbed = obed.readlines()	
	newbed = open(os.path.join(os.path.dirname(bed), '_gene_mappings.'.join([os.path.basename(bed).replace('.'+fileformat,''),fileformat])),'w')
	newbed2 = open(os.path.join(os.path.dirname(bed), '_isoform_mappings.'.join([os.path.basename(bed).replace('.'+fileformat,''),fileformat])),'w')
	k = 0
	m = 0
	n = 0
	o = 0
	lrbed = len(rbed)
	bar = progressbar.ProgressBar(maxval=lrbed, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in rbed:
		#print([line])
		isoforms = []
		liney = line.strip().split('\t')
		newbed.write('\t'.join(liney[0:9]) + '\t' + '1' + '\t' + str(int(liney[2]) - int(liney[1])) + ',' + '\t' + '0,' + '\n')
		features = liney[12:]
		nfeatures = len(liney[12:])
		#print(features)
		for i in features:
			if not i == '':
				iso = i.split(':')
				size = iso[0].split('-')[1]
				sumsize = sum([int(x) for x in size.strip(',').split(',')])
				nfeatures = len([int(x) for x in size.strip(',').split(',')])
				newbed2.write('\t'.join(liney[0:2]) + '\t' + str(int(liney[1])+sumsize) +'\t'+ '\t'.join(liney[3:7]) + '\t' + str(int(liney[1])+sumsize) + '\t' + \
							  liney[8] + '\t' + str(nfeatures) + '\t' + size + '\t' + iso[1].strip() + '\n')
		k += 1
		bar.update(k)
	bar.finish()
	newbed.close()
	newbed2.close()
	
def pick_genes(file, fasta): #Use this for speed
	#This script takes a fasta file and a file with the names of reads u want to take from the fasta. It will create a dictionary for the fasta seqs and search read names against the fasta name. 
	#If a name in the file is found in the fasta as well, the corresponding sequence and the name are returned. Pay attention to the \n and > and fasta format, it has to be single line for seqs
	print('Starting to pick gene fastas')
	openfasta = open(fasta, 'r')
	lopenfasta = openfasta.readlines()
	openfile=open(file,'r')
	lopenfile = openfile.readlines()
	lopenfile = lopenfile[1:]
	l = len(lopenfile)
	path = os.path.dirname(file)
	name = os.path.basename(file)
	outa=open(path+'/'+name.replace('.bed','')+'.fasta','w')
	outb=open(path+'/'+name.replace('.bed','')+'_unfound','w')

	n = 0
	x = 0
	g = 1
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for gene in lopenfile:
		ref = '>'+gene.split('\t')[0]+'\n'
		start = int(gene.split('\t')[1])
		end = int(gene.split('\t')[2])
		if start == 0:
			g = 0
		if ref in lopenfasta:
			ind = lopenfasta.index(ref)
			outa.write(ref.strip()+':'+str(start)+'_'+str(end)+'\t'+gene.split('\t')[3]+'\n'+lopenfasta[ind+1][start-g:end]+'\n')
			#outa.write(ref.strip()+':'+str(start)+'_'+str(end)+'\t'+gene.split('\t')[3]+'\n'+lopenfasta[ind+1][start-g:end]+'\n')
			#lopenfasta.remove(lopenfasta[ind])
			#lopenfasta.remove(lopenfasta[ind])
			x += 1
		else:
			outb.write(ref)
		n+=1
		bar.update(n)
	bar.finish()
	print('Discovered ' + str(x) + ' sequences' + ' which is %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*x/l),1))))
	
dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/capitata/rna_seq/Cc_E_5H/combined/minimap2_all/alignQC.ouput2/data_extract'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/SQANTI_tofu_min2_polyAtrim/novel_genes_curated'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/tofu/unaligned/unique_10X_genes/10x_only_genes/alignqc/gmap'

#collapse_genes(dr+'/'+'best.sorted_edited_resorted_over2689817.bed','bed')
#collapse_genes(dr+'/'+'best.sorted_edited_resorted_1008146-1624842.bed','bed')
#collapse_genes(dr+'/'+'best.sorted_edited_resorted_over705080.bed','bed')
#collapse_genes(dr+'/'+'best.sorted_edited_resorted_2144171-2689817.bed','bed')
#collapse_genes(dr+'/'+'best.sorted_edited_resorted_over2689817.bed','bed')
#collapse_genes(dr+'/'+'best.sorted_edited_resorted_over1700354.bed','bed')

#pick_genes(dr+'/'+'best.sorted2_sorted_1-308219_collapsed_genes_2.bed',dr+'/'+'olivefly_supernova_v1p1_corrvalid_group_gt10_500000bc_74X_seed0.pseudohap.1_ed_sl.fasta')
#pick_genes(dr+'/'+'best.sorted2_sorted_308220-614351_collapsed_genes_2.bed',dr+'/'+'olivefly_supernova_v1p1_corrvalid_group_gt10_500000bc_74X_seed0.pseudohap.1_ed_sl.fasta')
pick_genes(dr+'/'+'best.sorted.bed',dr+'/'+'olivefly_supernova_v1p1_corrvalid_group_gt10_500000bc_74X_seed0.pseudohap.1_sl2.fasta')

#parse_bed(dr+'/'+'best.sorted.bed','bed')

#os.system('cd %s && cat best.sorted_edited_resorted_1-519205_collapsed_genes.bed best.sorted_edited_resorted_519206-1002788_collapsed_genes.bed best.sorted_edited_resorted_over1002788_collapsed_genes.bed > best.sorted_edited_resorted_collapsed_genes.bed' %(dr))
#gene_mapping(dr+'/'+'best.sorted_edited_resorted_collapsed_genes.bed', 'bed')

#collapse_genes(dr+'/'+'best.sorted_resorted.bed','bed')
#collapse_genes(dr+'/'+'best.sorted_resorted_collapsed_genes.bed','bed')
