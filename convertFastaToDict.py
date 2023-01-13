#!/usr/bin/env python3

import re, os, sys, shelve #, patternCount5
import progressbar

def fasta_to_dict(fasta):
	#This function takes a fasta file with single line sequence and returns a dictionary with sequence name and corresponding read
	fastaNameSeq = {}
	fastanames = []
	fastaseqs = []
	openfasta = open(fasta, 'r')
	lopenfasta = openfasta.readlines()
	
	for i in lopenfasta:
		if re.match(r'^>', i):
			fastanames = fastanames + [i[0:len(i)-1]]
		elif re.match(r'^\w+', i):
			fastaseqs = fastaseqs + [i[0:len(i)-1]]

	for i in range(len(fastanames)):
		fastaNameSeq.setdefault(fastanames[i], fastaseqs[i])
			
	#print(fastaNameSeq)
	return(fastaNameSeq)
	
def contig(qry, ref):
	
	qryNameSeq = {}
	qrynames = []
	qryseqs = []
	qryfasta = open(qry, 'r')
	lqryfasta = qryfasta.readlines()
	
	refNameSeq = {}
	refnames = []
	refseqs = []
	reffasta = open(ref, 'r')
	lreffasta = reffasta.readlines()
	
	for i in lqryfasta:
		if re.match(r'^>\w+\d$', i):
			qrynames = qrynames + [i[0:len(i)-1]]
		elif re.match(r'^\w+', i):
			qryseqs = qryseqs + [i[0:len(i)-1]]
			
	for i in range(len(qrynames)):
		qryNameSeq.setdefault(qrynames[i], qryseqs[i])
		
	for j in lreffasta:
		if re.match(r'^>\w+\d$', j):
			#contigNameSeq.setdefault(i, '')
			refnames = refnames + [j[0:len(j)-1]]
		elif re.match(r'^\w+', j):
			refseqs = refseqs + [j[0:len(j)-1]]
	for j in range(len(refnames)):
		refNameSeq.setdefault(refnames[j], refseqs[j])
		
	refqrydict = shelve.open('/home/abayega/Desktop/hmm/refqrydict.txt')
	refqrydict['refdict'] = refNameSeq
	refqrydict['qrydict'] = qryNameSeq
	refqrydict.close()
	
	reffasta.close()
	qryfasta.close()
	
def openFiles():
	refdict = shelve.open('/gs/project/wst-164-ab/anthony/pacBio_data/olive_fly/jsa/assemblyAnalysis/refqrydict.txt')
	qrydict = shelve.open('/gs/project/wst-164-ab/anthony/pacBio_data/olive_fly/jsa/assemblyAnalysis/refqrydict.txt')
	
	#patternCount5.patternCount(qrydict['qrydict'], refdict['refdict'])
	patternCount(qrydict['qrydict'], refdict['refdict'])
	
def patternCount(qry, ref):
	#This code takes 2 dictionaries each containing keys as sequence names and values as sequences. It checks for unique
	#sequences in ref not in qry and prints em out as fasta, then it checks unique sequences in qry not in ref and prints em out in fasta
	for i,j in ref.items():
		refnames = ref.keys()
		refseq = ref.values()
	for i,j in qry.items():
		qrynames = qry.keys()
		qryseq = qry.values()
		
	nrefnames = []
	nrefseq = []
	nqrynames = []
	nqryseq = []
	
	summary1 = open('ref_summary.txt', 'a')
	summary2 = open('qry_summary.txt', 'a')
	
	for i in range(len(refseq)):
		if refseq[i] not in qryseq:
			#remove sequence and sequence name from lists; ref and qry
			nrefseq = nrefseq + [refseq[i]]
			nrefnames = nrefnames + [refnames[i]]
	for i in range(len(nrefnames)):
		summary1.write('>{0}'.format(nrefnames[i]))
		summary1.write('\n')
		summary1.write(nrefseq[i])
		summary1.write('\n')
		summary1.close()
		
	for i in range(len(qryseq)):
		if qryseq[i] not in refseq:
			#remove sequence and sequence name from lists; ref and qry
			nqrynames = nqrynames + [qrynames[i]]
			nqryseq = nqryseq + [qryseq[i]]
	for i in range(len(nqryseq)):
		summary2.write('>{0}'.format(nqrynames[i]))
		summary2.write('\n')
		summary2.write(nqryseq[i])
		summary2.write('\n')
		summary2.close()
		
def patternCount4(qry, ref):
	
	refseq = ref.values()
	summary = open('summary.txt', 'a')
	n = 0 #for looping thru qry seqs
	m = 0 #for looping thru ref seqs
	
	while m < (len(refseq)):
		for i,j in qry.items():
			if len(j) < len(refseq[m]):
				iterations = len(refseq[m]) - len(j) + 1
				for x in range(iterations):
					slide = str(refseq[m][n:n+len(j)])
					if slide == j:
						pos = n + 1
						writ = [ref.keys()[ref.values().index(refseq[m])]] + [len(refseq[m])] + [pos] + [i] + [len(j)]
						summary.write(str(writ))
						summary.write('\n')
					n = n + 1
				n = 0
			m = m + 1
	summary.close()
			
#ref = ['GACCATACTGA', 'GACCATACTGGACCATACTG' 'GACCGACCAGACCAGACCATGACCATACTGGACCATACTGGACCATACTGGACCATACTGACTGGACCATACTGTACTGTACGACCATACTGTGATACTGA',]
#qry = ['CATGACCATACTG', 'GACCATACTG']

#patternCount(qry, ref)

#fasta_to_dict('/mnt/parallel_scratch_mp2_wipe_on_december_2017/ioannisr/banthony/analysis/illumina/olivefly/rna_seq/hisat2_alignments/male_to_female_map/2.fasta')
look_up_given_reads('/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/ercc/gmap/mapped_reads',\
					'/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/ercc/gmap/Bo.E.5H_all_rnaseq_combined_pass.trimmed_stranded2.fasta')
#look_up_given_reads_lengths('sequences.mpa-format.labels.sorted4','/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/illumina/olivefly/kraken_male_genome/14.txt')