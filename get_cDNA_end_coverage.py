#!/usr/bin/env python3

import re, os, progressbar

def get_chromosome_features(bedtool_out, bed):
	newfile = open(bedtool_out+'_gene_coverage', 'w')
	newfile.write('gene\t'+ 'len\t' + 'cov_5\t' + 'cov_gene\t' + 'cov_3\t' + 'rat_cov_5\t' + 'rat_cov_3\n')
	dirc = os.path.dirname(bedtool_out)	
	
	rbed = open(bed)
	lbed = rbed.readlines()
	
	k = 0
	n = 0
	m = 0
	bar = progressbar.ProgressBar(maxval=len(lbed), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for gene in lbed:
		gene_name = gene.split('\t')[0]
		gene_len = int(gene.split('\t')[2].strip())
		os.system("cd %s && grep '%s' %s > gene.delete.bed" %(dirc,gene_name,bedtool_out))
		rgbed = open(dirc+'/'+'gene.delete.bed')
		lgbed = rgbed.readlines()
		med_5_cov = int(lgbed[49].split('\t')[4].strip())
		med_g_cov = int(lgbed[int("{0:.0f}".format(round(gene_len/2,0)))].split('\t')[4].strip())
		med_3_cov = int(lgbed[gene_len-49].split('\t')[4].strip())
		if med_g_cov != 0:
			newfile.write(gene_name + '\t' + str(gene_len) + '\t' + str(med_5_cov) + '\t' + str(med_g_cov) + '\t' + str(med_3_cov) + '\t' + str(round(med_5_cov/med_g_cov,2)) + '\t' + str(round(med_3_cov/med_g_cov,2)) + '\n')
		k += 1
		bar.update(k)
	bar.finish()
	newfile.close()

dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/Bo_E_1H/C010_10_171117/minimap2_ncbi_gff_all_reads/bedtools'
get_chromosome_features(dr+'/'+'positive_strand_genes_common_positive_Bo_E_1H_C010_10_bedtools_coverage.out.bed', dr+'/'+'positive_strand_genes_common_positive.bed')
		