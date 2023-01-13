#!/usr/bin/env python3

import os, re
import progressbar

def write_gene_list(gtf):
	a=b=c=d=e=f=g = 0
	genes=[]
	dirc = os.path.dirname(gtf)
	nfile1 = open(dirc+'/genelist', 'w')
	for line in open(gtf):
		a=line.strip().split('\t')
		if len(a)>6:
			if a[2]=='exon':
				start=int(a[3])
				end=int(a[4])
				chromosome=a[0]
				try:
					gene_name=a[8].split('gene_id "')[1].split('"')[0]+'_'+chromosome					
					if gene_name not in genes:
						genes+=[gene_name]
						nfile1.write(gene_name.strip()+'\n')
						c+=1
				except:
					gene_name=a[8].split('transcript_id "')[1].split('"')[0]+'_'+chromosome
					if gene_name not in genes:
						genes+=[gene_name]
						nfile1.write(gene_name.strip()+'\n')
						d+=1
					#print(line)
					#sys.exit(1)
	nfile1.close()
	print('\nWritten '+str(c)+' genes and '+str(d)+' transcript IDs\n'+'\nTotal genes and transcript IDs is '+str(len(genes))+'\n')
					
def count_chromosome_genes(gtf, psl):	
	o = 0
	k = 0
	badpsl = [] 
	'''
	print('\nChecking alignment file for errors, hang on a sec...\n')
	lenpsl = int(os.popen('wc -l %s' %(psl)).read().split()[0])
	bar = progressbar.ProgressBar(maxval=lenpsl, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for i in open(psl):
		if len(i.split('\t')) > 30 or len(i.split('\t')) < 20 or len(i.split('\t')[20].split(',')) == 0:
			badpsl += [i]			
			o += 1
		k += 1
		bar.update(k)
	bar.finish()
	
	if o >= 1:
		print('\nThere seems to be ' + str(o) + ' number of bad sequences. These are writen in edited version of psl file...\n\n')
		newpsl = open(os.path.dirname(psl)+'/'+os.path.basename(psl)+'_edited', 'w')
		for x in badpsl:
			newpsl.write(x)
		newpsl.close()
	else:
		print('\nThe alignment file should be fine, moving on...\n\n')
	'''
		
	dirc = os.path.dirname(gtf)
	if not os.path.exists(dirc+'/output2'):
		os.makedirs(dirc+'/output2')
	print('Counter should be equal to ' + str(os.popen('wc -l %s' %(os.path.dirname(gtf)+'/2D_trimmed_l_gmapoutput_filtered.psl')).read().split()[0]))
	os.system("cd %s && awk '//{print($1)}' %s | sort -u > chromosome.txt" %(os.path.dirname(gtf),gtf)) #U can supply your own file with chromosomes to check
	print('Done generating chromosomes')
	chromosomes = int(os.popen('wc -l %s' %(os.path.dirname(gtf)+'/chromosome.txt')).read().split()[0])
	k = 1
	for line in open(os.path.dirname(gtf)+'/chromosome.txt'):
		os.system('cd %s && rm new_file.gtf' %(os.path.dirname(gtf)))
		os.system('cd %s && rm new_alignments.psl' %(os.path.dirname(gtf)))
		print('Chromosome ' + str(k) + ' of ' + str(chromosomes))
		#create chromosome gtf
		os.system('cd %s && grep %s %s >> new_file.gtf' %(os.path.dirname(gtf), line.strip(), gtf))
		print('Done generating chromosome ' + line.strip() + ' gtf')
		#get chromosome alignments
		os.system('cd %s && grep %s %s >> new_alignments.psl' %(os.path.dirname(gtf), line.strip(), psl))
		print('Done generating chromosome ' + line.strip() + ' alignments')
		#get gene list
		newfile = open(os.path.dirname(gtf)+'/chromosome_genelist.txt', 'w')
		for linei in open(os.path.dirname(gtf)+'/new_file.gtf'): #I changed the code here coz previously it was "for linei in open(gtf):" 
			linex = linei.split()
			if 'gene_id' in linex:
				newfile.write(linex[linex.index('gene_id')+1].strip(';"') + '\n')
		newfile.close()
		os.system("cd %s && sort -u chromosome_genelist.txt > chromosome_genelist2.txt" %(os.path.dirname(gtf)))
		print('Done generating chromosome ' + line.strip() + ' gene list')
		os.system('cd %s && bash mandalorion.sh' %(os.path.dirname(gtf)))
		os.system('cat %s >> %s' %(dirc+'/'+'output/Gene_Expression.txt',dirc+'/'+'output2/Gene_Expression_all.txt' ))
		print('Done generating chromosome ' + line.strip() + ' Mandalorion gene counts. Moving on...')
		k += 1
		#if k == 70:
		#	break


def get_chromosome_features(gtf, psl):
	dirc = os.path.dirname(gtf)
	if not os.path.exists(dirc+'/output2'):
		os.makedirs(dirc+'/output2')
	print('Counter should be equal to ' + str(os.popen('wc -l %s' %(os.path.dirname(gtf)+'/2D_trimmed_l_gmapoutput_filtered.psl')).read().split()[0]))
	os.system("cd %s && awk '//{print($1)}' %s | sort -u > chromosome.txt" %(os.path.dirname(gtf),gtf))
	print('Done generating chromosomes')
	chromosomes = int(os.popen('wc -l %s' %(os.path.dirname(gtf)+'/chromosome.txt')).read().split()[0])
	k = 1
	for line in open(os.path.dirname(gtf)+'/chromosome.txt'):
		os.system('cd %s && rm -r parsed_reads' %(dirc))
		os.system('cd %s && rm new_file.gtf' %(os.path.dirname(gtf)))
		os.system('cd %s && rm new_alignments.psl' %(os.path.dirname(gtf)))
		print('Chromosome ' + str(k) + ' of ' + str(chromosomes))
		#create chromosome gtf
		os.system('cd %s && grep %s %s >> new_file.gtf' %(os.path.dirname(gtf), line.strip(), gtf))
		print('Done generating chromosome ' + line.strip() + ' gtf')
		#get chromosome alignments
		os.system('cd %s && grep %s %s >> new_alignments.psl' %(os.path.dirname(gtf), line.strip(), psl))
		print('Done generating chromosome ' + line.strip() + ' alignments')
		#get gene list
		newfile = open(os.path.dirname(gtf)+'/chromosome_genelist.txt', 'w')
		for linei in open(os.path.dirname(gtf)+'/new_file.gtf'): #I changed the code here coz previously it was "for linei in open(gtf):"
			linex = linei.split()
			if 'gene_id' in linex:
				newfile.write(linex[linex.index('gene_id')+1].strip(';"') + '\n')
		newfile.close()
		os.system("cd %s && sort -u chromosome_genelist.txt > chromosome_genelist2.txt" %(os.path.dirname(gtf)))
		print('Done generating chromosome ' + line.strip() + ' gene list')
		os.system('cd %s && bash mandalorion.sh' %(os.path.dirname(gtf)))
		for file in os.listdir(dirc+'/output'):
			os.system('cat %s >> %s' %(dirc+'/output/'+file, dirc+'/output2/'+file))
		print('Done generating chromosome ' + line.strip() + ' Mandalorion gene counts. Moving on...')
		k += 1	

dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/Bo_E_2H/C010_09_141117/mandalorion_MoY'
dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/Bo_E_3H/C010_08_071117/mandalorion_MoY'
dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/Bo_E_4H/C010_06_191017/mandalorion_MoY'
dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/Bo_E_5H/combined/mandalorion_MoY'
dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/Bo_E_6H/C010_07_4_041117/mandalorion_MoY'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_1H/C010_10_171117/mandalorion'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_2H/C010_09_141117/mandalorion'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_3H/C010_08_071117/mandalorion'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_4H/C010_06_191017/mandalorion'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_5H/combined/mandalorion'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_6H/C010_07_4_041117/mandalorion'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_Heads_Female/C010_12_01_120518/mandalorion'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_Heads_Male/C010_12_02_140518/mandalorion'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/capitata/rna_seq/Cc_E_5H/C010_03C_301017/mandalorion'

#count_chromosome_genes(dr+'/'+'100279_2148_2850.gtf', dr+'/'+'2D_trimmed_l_gmapoutput_filtered.psl')
#get_chromosome_features(dr+'/ERCC92_ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_gffread_sorted_edited_over37.gtf', dr+'/2D_trimmed_l_gmapoutput_filtered.psl')
write_gene_list(dr+'/'+'ERCC92_GCF_000347755.3_Ccap_2.1_genomic.gtf')

