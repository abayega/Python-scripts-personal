#!/usr/bin/env python3

from __future__ import division
import re, os, sys, shelve, progressbar
from baseconv import base2
import statistics

def checkfiles(qry):
	newfile = open('newfile.txt', 'w')
	newfile2 = open('newfilex.txt', 'w')
	writ = 'read_identity' + '\t'+ 'Contig_aligned' + '\t' + 'read_length' + '\t'+ 'ref_length' + '\t'+ 'align_identity' + \
	'\t'+ 'read_identity' + '\t' + 'insertion' + '\t' + 'deletion' + '\t' + 'identical-bases' + '\t' + 'mismatches' + '\t' + '\n'
	writ = 'read_identity' + '\t' + 'SAM-FLAG' + '\t' + 'ref-ID' + '\n'
	
	newfile.write(str(writ))
	samfile = open(qry, 'r')
	samline = samfile.readlines()
	for i in samline:
		if not i.startswith('@'):
			sam_analyse(i)

def purse_cigar(qry):
	#This script parses sam file lines and returns several attributes about the CIGAR
	print('Starting purse_cigar analysis')
	newfile = open(os.path.join(os.path.dirname(qry), '_'.join([os.path.basename(qry).replace('.sam',''),'with_alignments.txt'])), 'w')
	newfile2 = open(os.path.join(os.path.dirname(qry), '_'.join([os.path.basename(qry).replace('.sam',''),'without_alignments.txt'])), 'w')
	writ = 'read_name\t' + 'Ref_aligned\t' + 'read_length\t' + 'alignment_length\t' + 'aligned_length\t' + 'insertions\t' + \
	'deletions\t' + 'indels\t' + 'median_ins\t' + 'median_del\t' + 'max_ins\t' + 'max_del\t' + 'read_aligned\t' + 'ref_length\t' + 'readln_check' + '\n'
	writ2 = 'read_identity\t' + 'SAM_FLAG\t'+ 'ref-ID' + '\n'
	newfile.write(str(writ))
	newfile2.write(str(writ2))
	
	#samfile = open(qry, 'r')
	#samline = samfile.readlines()
	samline = int(os.popen('wc -l %s' %(qry)).read().split()[0])
	badlines = 0
	rdsaln = 0
	noalign = 0
	p=0
	bar = progressbar.ProgressBar(maxval=len(samline), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	#for i in samline:
	for line in open(qry):
		alignment = i.split('\t')
		try:
			extra_fields = alignment[11:]
			readln1 = 'unknown'
			refln = 'unknown'
			for i in extra_fields:
				if 'ZQ:i:' in i:
					readln1 = int(i.lstrip('ZQ:i:'))
				elif 'ZR:i:' in i:
					refln = int(i.lstrip('ZR:i:'))
					
			readln2 = len(alignment[9])
			
			#cigar length
			pat = re.compile(r'(\D)', re.I)
			npat = pat.sub(r'+', alignment[5])
			npatlist = npat.strip('+').split('+')
			sumnpat = 0
			for x in npatlist:
				sumnpat += int(x)
			
			#soft clipping
			sumsoft = 0
			spat = re.compile(r'(\d+)S', re.I)
			soft = spat.findall(alignment[5])
			for x in soft:
				sumsoft += int(x)
				
			#matches (identical + mismatches)
			summatches = 0
			matpat = re.compile(r'(\d+)M', re.I)
			mpat = matpat.findall(alignment[5])
			for x in mpat:
				summatches += int(x)
			
			#insertions
			sumins = 0
			ipat = re.compile(r'\w(\d+)I', re.I)
			ins = ipat.findall(alignment[5])
			ins2 = []
			for x in ins:
				ins2+=[int(x)]				
			sumins = sum(ins2)
			medins = statistics.median(ins2)
			maxins = max(ins2)
			
			#deletions
			sumdeln = 0
			dpat = re.compile(r'\w(\d+)D', re.I)
			deln = dpat.findall(alignment[5])
			deln2 = []
			for x in deln:
				deln2+=[int(x)]				
			sumdeln = sum(deln2)
			meddeln = statistics.median(deln2)
			maxdeln = max(deln2)
			
			align_len = sumnpat - sumsoft - sumdeln  #The portion of the raw query read alignment
			
			#alignment and read identity
			readln = sumnpat - sumdeln
			if readln == readln1:
				check = 1
			elif readln == readln2:
				check = 11
			elif readln1 == readln2:
				check = 111
			else:
				check = 0
			alignmentln = sumnpat - sumsoft
			alignedln = alignmentln - sumdeln
			sumerrors = alignmentln # alignmentln == sumerrors == alignment_string_length
			alignment_string_length = sumerrors	
			
			read_ins = 100 * sumins/sumerrors
			read_deln = 100 * sumdeln/sumerrors
			indels = 100 * (sumins + sumdeln)/sumerrors
			perc_aligned = 100 * align_len/readln			

			#write results
			writ = alignment[0] + '\t' + alignment[2] + '\t' + str(readln) + '\t' + str(alignmentln) + '\t' + str(alignedln) + '\t' + str(round(read_ins,1)) + '\t' + str(round(read_deln,1)) + \
			'\t' + str(round(indels,1)) + '\t' + str(medins) + '\t' + str(meddeln) + '\t' + str(maxins) + '\t' + str(maxdeln) + '\t' + str(round(perc_aligned,0)) + '\t' + str(refln) + '\t' + str(check) + '\n'
			newfile.write(str(writ))
			rdsaln += 1
			
		except:
			if not alignment[0].startswith('@'):
				writ2 = alignment[0] + '\t'+ alignment[1] + '\t'+ alignment[2] + '\n'
				newfile2.write(str(writ2))
				noalign += 1
			else:
				badlines += 1
		p+=1
		bar.update(p)
	bar.finish()
	print('There were ' + str(badlines) + ' lines starting with @__ and ' + str(noalign) + ' lines with no alignments and are recorded.' + '\n' + 'Percentage of lines showing alignments = ' + str(round(rdsaln/(rdsaln + noalign) * 100,0)))
	print('Finised, bye')
	newfile.close()
	newfile2.close()
	
def sam_analyse_hisat2_mem(qry):
	#This script takes a multifasta .sam file output from bwa-mem and return read_name, 
	#ref_aligned, read_length, ref_length, align_identity, read_identity, insertion rate, 
	#deletion rate, identical-bases, mismatches, %_read-aligned, and sam-FLAG in a tab delimited format. 
	#Two files are returned, one with reads that have an alignment and another with reads that have no alignment.
	#newfile = open(os.path.join(os.path.dirname(qry), 'ecoli_1-15kb.contigs_with_alignments.txt'), 'w')
	#newfile2 = open(os.path.join(os.path.dirname(qry), 'ecoli_1-15kb.contigs_without_alignments.txt'), 'w')
	print('Starting bwa sam file analysis')
	newfile = open(os.path.join(os.path.dirname(qry), '_'.join([os.path.basename(qry).replace('.sam',''),'with_alignments.txt'])), 'w')
	newfile2 = open(os.path.join(os.path.dirname(qry), '_'.join([os.path.basename(qry).replace('.sam',''),'without_alignments.txt'])), 'w')
	writ = 'read_name' + '\t' + 'Ref_aligned' + '\t' + 'read_length' + '\t' + 'align_identity' + '\t' + 'read_identity' + '\t' + 'insertion' + \
	'\t' + 'deletion' + '\t' + 'identical_bases' + '\t' + 'mismatches' + '\t' + 'perc_read_aligned' + '\t' + 'SAM_FLAG' + '\t' + 'indels' + '\t' + \
	'aligned_identity' + '\t' + 'sumerrors' + '\t' + 'aligned_string_length' + '\t' + 'A' + '\t' + 'C' + '\t' + 'G' + '\t' + 'T' + '\t' + 'align_score' + '\t' + 'flag' + '\n'
	writ2 = 'read_identity' + '\t' + 'SAM_FLAG' + '\t'+ 'ref_ID' + '\n'
	newfile.write(str(writ))
	newfile2.write(str(writ2))
	
	#with open(qry, encoding="latin-1") as samfile:
	samfile = open(qry, 'r')
	samline = samfile.readlines()
	#samlines = int(os.popen('wc -l %s' %(qry)).read().split()[0])
	badlines = 0
	rdsaln = 0
	noalign = 0
	p=0
	bar = progressbar.ProgressBar(maxval=len(samline), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for i in samline:
	#for i in open(qry):
		alignment = i.split('\t')
		try:
			extra_fields = alignment[11:]
			align_score = [i.strip('AS:i:') for i in extra_fields if 'AS:i:' in i]
			if len(align_score) > 0:
				align_score = align_score[0]
			else:
				align_score = ''
			
			if alignment[1].isdecimal():
				sam_flag = sam_flag_debug(alignment[1])
			else:
				sam_flag = ''
				
			#if alignment[1].is and alignment[1] != 4 and alignment[2] != '*':
			pat = re.compile(r'(\D)', re.I)
			npat = pat.sub(r'+', alignment[5])
			npatlist = npat.strip('+').split('+')
			sumnpat = 0
			for x in npatlist:
				sumnpat += int(x)
			
			#matches (identical + mismatches)
			summatches = 0
			matpat = re.compile(r'(\d+)M', re.I)
			mpat = matpat.findall(alignment[5])
			for x in mpat:
				summatches += int(x)
			
			#Ns gapped alignments
			sumNs = 0
			npat = re.compile(r'(\d+)N', re.I)
			Ns = npat.findall(alignment[5])
			for x in Ns:
				sumNs += int(x)
			
			#soft clipping
			sumsoft = 0
			spat = re.compile(r'(\d+)S', re.I)
			soft = spat.findall(alignment[5])
			for x in soft:
				sumsoft += int(x)
			
			#insertions
			sumins = 0
			ipat = re.compile(r'\w(\d+)I', re.I)
			ins = ipat.findall(alignment[5])
			for x in ins:
				sumins += int(x)
			
			#deletions
			sumdeln = 0
			dpat = re.compile(r'\w(\d+)D', re.I)
			deln = dpat.findall(alignment[5])
			for x in deln:
				sumdeln += int(x)
			align_len = sumnpat - sumsoft - sumdeln  
			align_len = summatches + sumins #The portion of the raw query read alignment
			
			#sequence matches from MD
			MD = [i.strip('MD:Z:') for i in extra_fields if 'MD:Z:' in i]
			MD = MD[0]
			if len(MD) > 0 and not MD.isdecimal():				
				A = MD.count("A")
				T = MD.count("T")
				C = MD.count("C")
				G = MD.count("G")				    
				mismatches = A + C + G + T #These are reference bases which are wrongly called in the querry, meaning the querry was another base
			else:
				A = 0
				T = 0
				C = 0
				G = 0
				mismatches = 0
			
			align_seqmatch = summatches - mismatches			
			
			#alignment identity
			sumerrors = sumins + sumdeln + align_seqmatch +  mismatches			
			readln = sumnpat - sumdeln - sumNs
			
			alignment_string_length = summatches + sumins + sumdeln
			
			align_identity = 100 * align_seqmatch/alignment_string_length
			read_identity = 100 * align_seqmatch/readln # = overall_identity
			aligned_identity = 100 * (align_seqmatch / (summatches + sumins))			
			
			read_ins = 100 * sumins/sumerrors
			read_deln = 100 * sumdeln/sumerrors
			read_seqmat = 100 * align_seqmatch/sumerrors
			read_mismat = 100 * mismatches/sumerrors
			perc_aligned = 100 * align_len/readln
			indels = 100 * (sumins + sumdeln)/sumerrors
			
			if mismatches == 0:
				mismatches = 1

			#write results
			writ = alignment[0] + '\t' + alignment[2] + '\t' + str(readln) + '\t' + str("{0:.0f}".format(round(align_identity,0))) + '\t' + \
			str("{0:.0f}".format(round(read_identity,0))) + '\t' + str("{0:.1f}".format(round(read_ins,1))) + '\t' + str("{0:.1f}".format(round(read_deln,1))) + \
			'\t' + str("{0:.1f}".format(round(read_seqmat,1))) + '\t' + str("{0:.1f}".format(round(read_mismat,1))) +'\t' + str("{0:.0f}".format(round(perc_aligned,0))) + \
			'\t' + alignment[1] + '\t' + str("{0:.1f}".format(round(indels,1))) + '\t' + str("{0:.1f}".format(round(aligned_identity,1))) + '\t' + \
			str("{0:.0f}".format(round(sumerrors,0))) + '\t' + str("{0:.0f}".format(round(alignment_string_length,0))) + '\t' + \
			str("{0:.0f}".format(round(100*A/mismatches,0))) + '\t' + str("{0:.0f}".format(round(100*C/mismatches,0))) + '\t' + \
			str("{0:.0f}".format(round(100*G/mismatches,0))) + '\t' + str("{0:.0f}".format(round(100*T/mismatches,0))) + '\t' + align_score +  '\t' + sam_flag + '\n'            
			#writ = alignment[0] + '\t'+ alignment[2] + '\t'+ str(align_identity) + '\t'+ str(read_identity) + '\t'+ str(readln) + '\t'+ str(refln) + '\t' + str(read_ins) + '\t' + str(read_deln) + '\t' + alignment[1] + '\n'
			#if not sam_flag == '256' and not sam_flag == '4':
			if sam_flag != '256' and sam_flag != '4':
				newfile.write(str(writ))
			#newfile.write(str(writ))
			rdsaln += 1
			#print(rdsaln)
			
		except:
			if not alignment[0].startswith('@'): #and alignment[1] != 4 and alignment[2] != '*':
				writ2 = alignment[0] + '\t'+ alignment[1] + '\t'+ alignment[2] + '\n'
				newfile2.write(str(writ2))
				noalign += 1
			else:
				badlines += 1
				#continue
		p+=1
		bar.update(p)
	bar.finish()
	#print('There were ' + str(badlines) + ' reads without alignments in the file.' + '\n' + 'Percentage of reads with alignments = ' + str("{0:.0f}".format(round(rdsaln/(badlines + rdsaln) * 100,0))))
	print('There were ' + str(badlines) + ' lines starting with @__ and ' + str(noalign) + ' lines with no alignments and are recorded.' + '\n' + 'Percentage of lines showing alignments = ' + str("{0:.0f}".format(round(rdsaln/(badlines + rdsaln + noalign) * 100,0))))
	print('Finised, bye')
	newfile.close()
	newfile2.close()
	
def sam_analyse_bwa_mem2(qry):
	#This script takes a multifasta .sam file output from bwa mem and return read_name, 
	#ref_aligned, read_length, ref_length, align_identity, read_identity, insertion rate, 
	#deletion rate, identical-bases, mismatches, %_read-aligned, and sam-FLAG in a tab delimited format. 
	#Two files are returned, one with reads that have an alignment and another with reads that have no alignment.
	#newfile = open(os.path.join(os.path.dirname(qry), 'ecoli_1-15kb.contigs_with_alignments.txt'), 'w')
	#newfile2 = open(os.path.join(os.path.dirname(qry), 'ecoli_1-15kb.contigs_without_alignments.txt'), 'w')
	print('Starting bwa sam file analysis')
	newfile = open(os.path.join(os.path.dirname(qry), '_'.join([os.path.basename(qry).replace('.sam',''),'with_alignments1.txt'])), 'w')
	newfile2 = open(os.path.join(os.path.dirname(qry), '_'.join([os.path.basename(qry).replace('.sam',''),'without_alignments1.txt'])), 'w')
	writ = 'read_name' + '\t' + 'SAM_FLAG' + '\t' + 'ref_name' + '\t' + 'map_quality' + '\t' + 'gapped_map' + '\t' + 'unique_map' + '\t' + 'multi_map' + \
	 '\t' + 'perfect_match' + '\t' + 'mismatches' + '\t' + 'read_seq'+ '\t'  + 'sam_flags' + '\n'
	writ2 = 'read_identity' + '\t' + 'SAM_FLAG' + '\t'+ 'ref-ID' + '\n'
	newfile.write(str(writ))
	newfile2.write(str(writ2))
	
	samfile = open(qry, 'r')
	samline = samfile.readlines()
	badlines = 0
	rdsaln = 0
	noalign = 0
	p=0
	bar = progressbar.ProgressBar(maxval=len(samline), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for i in samline:
		alignment = i.split('\t')
		try:
			extra_fields = alignment[11:]
			MD = [i.strip('MD:Z:') for i in extra_fields if 'MD:Z:' in i]
			if len(MD) > 0:
				MD = MD[0]
			else:
				MD = ''
			if 'XT:A:U' in extra_fields:
				unique_map = 'yes'
			else:
				unique_map = 'unknown'
			if 'XT:A:R' in extra_fields:
				multi_map = 'yes1'
			else:
				multi_map = 'unknown1'			
			if 'NM:i:0' in extra_fields or 'XM:i:0' in extra_fields:
				perfect_match = 'yes2'
				mismatches = 0
			elif MD.isdecimal():
				perfect_match = 'yes2'
				mismatches = 0
			else:
				perfect_match = 'no'
				mismatches = 'unknown3'
			for x in extra_fields:
				if 'NM:i:' in x:
					mismatches = x.strip('NM:i:')
				elif 'XM:i:' in x:
					mismatches = x.strip('XM:i:')
					
			if 'I' in alignment[5] or 'D' in alignment[5]:
				gapped_map = 'yes3'
			else:
				gapped_map = 'no'
			if alignment[1].isdecimal():
				sam_flag = sam_flag_debug(alignment[1]) 
			#write results			
			writ = alignment[0] + '\t' + alignment[1] + '\t' + alignment[2] + '\t' + alignment[4] + '\t' + gapped_map + '\t' + unique_map + '\t' + multi_map + '\t' + \
			perfect_match + '\t' + mismatches + '\t' + alignment[9] + '\t' + sam_flag + '\n'
			newfile.write(str(writ))
			rdsaln += 1
			#print(rdsaln)
			
		except:
			if not alignment[0].startswith('@'): #and alignment[1] != 4 and alignment[2] != '*':
				writ2 = alignment[0] + '\t'+ alignment[1] + '\t'+ alignment[2] + '\n'
				newfile2.write(str(writ2))
				noalign += 1
			else:
				badlines += 1
				#continue
		p+=1
		bar.update(p)
	bar.finish()
	#print('There were ' + str(badlines) + ' reads without alignments in the file.' + '\n' + 'Percentage of reads with alignments = ' + str("{0:.0f}".format(round(rdsaln/(badlines + rdsaln) * 100,0))))
	print('There were ' + str(badlines) + ' lines starting with @__ and ' + str(noalign) + ' unprocessed lines and are recorded.' + '\n' + 'Percentage of processed lines = ' + str("{0:.0f}".format(round(rdsaln/(badlines + rdsaln + noalign) * 100,0))))
	print('Finised, bye')
	newfile.close()
	newfile2.close()	
		
def sam_flag_debug(samflag):
	#This script takes a number or string of number from a SAM_FLAG and will return the debuged meanings of the number
	flag_dic = {0:'paired_read',1:'both_pairs_mapped',2:'unmapped',3:'mate_is_unmapped',40:'query_forward',41:'query_reverse',50:'mate_forward',51:'mate_reverse',\
				6:'first_sequence_of_a_pair',7:'second_sequence_of_a_pair',8:'not_primary_alignment',9:'sequence_or_alignment_does_not_pass_quality_controls',\
				10:'PCR_duplicate',11:'supp_align'} #8:'not_the_best_alignment_of_this_sequence'
	
	summary = []
	lis_encode_samflag = list(base2.encode(samflag)[::-1])
	lencode_samflag = len(lis_encode_samflag)
	#starting_code = len(lencode_samflag) - 1
	if lencode_samflag == 1 and lis_encode_samflag[0] == '0':
		return('query_forward_unpaired')

	for i in range(len(lis_encode_samflag)):		
		if i == 4 and int(lis_encode_samflag[i]) == 0:
			summary = summary + [flag_dic[int(i*10)]]
		elif i == 4 and int(lis_encode_samflag[i]) == 1:
			summary = summary + [flag_dic[int(i*10+1)]]
		elif i == 5 and int(lis_encode_samflag[i]) == 0:
			summary = summary + [flag_dic[int(i*10)]]
		elif i == 5 and int(lis_encode_samflag[i]) == 1:
			summary = summary + [flag_dic[int(i*10+1)]]
		elif i != 4 and i != 5 and int(lis_encode_samflag[i]) == 1:
			summary = summary + [flag_dic[int(i)]]			
	
	return(','.join(summary))

def sam_analyse_minimap2(qry): #, chimera):
	#This script takes a .sam file output from minimap2 and return read_name, 
	#ref_aligned, read_length, ref_length, align_identity, read_identity, insertion rate, 
	#deletion rate, identical-bases, mismatches, %_read-aligned, and sam-FLAG in a tab delimited format. 
	#Two files are returned, one with reads that have an alignment and another with reads that have no alignment.
	print('Starting minimap2 sam file analysis')
	newfile = open(os.path.join(os.path.dirname(qry), '_'.join([os.path.basename(qry).replace('.sam',''),'with_alignments.txt'])), 'w')
	newfile2 = open(os.path.join(os.path.dirname(qry),'_'.join([os.path.basename(qry).replace('.sam',''),'without_alignments.txt'])), 'w')
	#newfile3 = open(os.path.join(os.path.dirname(qry),'_'.join([os.path.basename(qry).replace('.sam',''),'good_alignments.sam'])), 'w')
	writ = 'read_name\t' + 'Ref_aligned\t' + 'read_length\t' + 'alignment_identity\t' + 'read_identity\t' + 'aligned_identity\t' + 'insertion\t' + \
	'deletion\t' + 'indels\t' + 'perc_read_aligned\t' + 'SAM-FLAG\t' + 'aligned_string_length\t' + 'chaining_score1\t' + 'chaining_score2\t' + 'aln_score\t' + \
	'ref_edit_distance\t' + 'strand_orientation\t' + 'aln_type\t' + 'readlength?\t' + 'flag\n'
	writ2 = 'read_identity\t' + 'SAM_FLAG\t' + 'ref-ID\n'
	newfile.write(str(writ))
	newfile2.write(str(writ2))
	
	samfile = open(qry, 'r')
	samline = samfile.readlines()
	#chimeric = open(chimera, 'r')
	#chimline = chimeric.readlines()
	#print('There are ' + str(len(chimline)) + ' chimeras')
	badlines = 0
	rdsaln = 0
	noalign = 0	
	p=0
	bar = progressbar.ProgressBar(maxval=len(samline), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in samline:		
		alignment = line.split('\t')
		try:
			if alignment[1].isdecimal():
				sam_flag = sam_flag_debug(alignment[1])
			else:
				sam_flag = ''
			
			extra_fields = alignment[11:]
			aln_type = '*'
			chaining_score1 = '*'
			chaining_score2 = '*'
			ref_edit_distance = '*'
			aln_score = '*'
			strand_orientation = '*'
			for i in extra_fields:
				if 'tp:A:' in i:
					aln_type = i.strip('tp:A:')
				elif 's1:i:' in i:
					chaining_score1 = i.strip('s1:i:')
				elif 's2:i:' in i:
					chaining_score2 = i.strip('s2:i:')
				elif 'NM:i:' in i:					
					ref_edit_distance = int(i.strip('NM:i:'))				
				elif 'AS:i:' in i:					
					aln_score = i.strip('AS:i:')
				elif 'ts:A:' in i:
					strand_orientation = i.strip('ts:A:')
				elif 'SA:Z:' in i and '+' in i:
					strand_orientation = '+'
				elif 'SA:Z:' in i and '-' in i:
					strand_orientation = '-'
				
			if alignment[1].isdecimal():
				sam_flag = sam_flag_debug(alignment[1])
			else:
				sam_flag = ''
			
			#if alignment[1].is and alignment[1] != 4 and alignment[2] != '*':
			pat = re.compile(r'(\D)', re.I)
			npat = pat.sub(r'+', alignment[5])
			npatlist = npat.strip('+').split('+')
			sumnpat = 0
			for x in npatlist:
				sumnpat += int(x)
			#align_len1 = sumnpat
			
			#matches (identical + mismatches)
			summatches = 0
			matpat = re.compile(r'(\d+)M', re.I)
			mpat = matpat.findall(alignment[5])
			for x in mpat:
				summatches += int(x)
			align_matches = summatches
				
			#soft clipping
			sumsoft = 0
			spat = re.compile(r'(\d+)S', re.I)
			soft = spat.findall(alignment[5])
			for x in soft:
				sumsoft += int(x)
			align_soft = sumsoft
			#align_len2 = sumnpat - align_soft
			
			#Ns gapped alignments
			sumNs = 0
			npat = re.compile(r'(\d+)N', re.I)
			Ns = npat.findall(alignment[5])
			for x in Ns:
				sumNs += int(x)
				
			#insertions
			sumins = 0
			ipat = re.compile(r'\w(\d+)I', re.I)
			ins = ipat.findall(alignment[5])
			for x in ins:
				sumins += int(x)
			align_ins = sumins
			
			#deletions
			sumdeln = 0
			dpat = re.compile(r'\w(\d+)D', re.I)
			deln = dpat.findall(alignment[5])
			for x in deln:
				sumdeln += int(x)
			align_deln = sumdeln
				
			align_len = align_matches  + align_ins			

			#alignment identity
			#readln = len(alignment[9])
			readln = sumnpat - sumNs - align_deln
			alignment_string_length = sumnpat - sumsoft - sumNs
			x = 1
			if sumnpat - align_deln - sumNs != readln:
				x = 0				
				
			alignment_identity = 100 - (100 * ref_edit_distance/alignment_string_length)
			read_identity = 100 - (100 * ref_edit_distance/readln)
			aligned_identity = 100 - (100 * ref_edit_distance/align_len)
			
			read_ins = 100 * align_ins/alignment_string_length
			read_deln = 100 * align_deln/alignment_string_length #This is what you should trust most together with percentage of read aligned
			indels = 100 * (align_ins + align_deln)/alignment_string_length
			perc_aligned = 100 * align_len/readln
			
			#write results
			#writ = alignment[0] + '\t' + alignment[2] + '\t' + str(readln) + '\t' + str("{0:.0f}".format(round(alignment_identity,0))) + '\t' + str("{0:.0f}".format(round(read_identity,0))) + '\t' + \
			#str("{0:.1f}".format(round(aligned_identity,1))) + '\t' + str("{0:.1f}".format(round(read_ins,1))) + '\t' + str("{0:.1f}".format(round(read_deln,1))) + '\t' + str("{0:.1f}".format(round(indels,1))) + '\t' + \
			#str("{0:.0f}".format(round(perc_aligned,0))) + '\t' + alignment[1] + '\t' + str("{0:.0f}".format(round(alignment_string_length,0))) + '\t' + \
			#chaining_score1 + '\t' + chaining_score2 + '\t' + aln_score + '\t' + str(ref_edit_distance) + '\t' + strand_orientation.strip('\n') + '\t' + aln_type + '\t' + str(x)  + '\t' +sam_flag + '\n'
			
			#newfile.write(str(writ))
			rdsaln += 1
			
			writ = alignment[0] + '\t' + alignment[2] + '\t' + str(readln) + '\t' + str("{0:.0f}".format(round(alignment_identity,0))) + '\t' + str("{0:.0f}".format(round(read_identity,0))) + '\t' + \
			str("{0:.1f}".format(round(aligned_identity,1))) + '\t' + str("{0:.1f}".format(round(read_ins,1))) + '\t' + str("{0:.1f}".format(round(read_deln,1))) + '\t' + str("{0:.1f}".format(round(indels,1))) + '\t' + \
			str("{0:.0f}".format(round(perc_aligned,0))) + '\t' + alignment[1] + '\t' + str("{0:.0f}".format(round(alignment_string_length,0))) + '\t' + \
			chaining_score1 + '\t' + chaining_score2 + '\t' + aln_score + '\t' + str(ref_edit_distance) + '\t' + strand_orientation.strip('\n') + '\t' + aln_type + '\t' + str(x)  + '\t' +sam_flag + '\n'
			newfile.write(str(writ))
			
			#if perc_aligned > 79 and perc_aligned < 101 and aln_type != 'S' and alignment[0]+'\n' not in chimline:
			#if perc_aligned > 50 and alignment_identity > 85 and aln_type != 'S':
			#	newfile3.write(line)
			#	writ = alignment[0] + '\t' + alignment[2] + '\t' + str(readln) + '\t' + str("{0:.0f}".format(round(alignment_identity,0))) + '\t' + str("{0:.0f}".format(round(read_identity,0))) + '\t' + \
			#	str("{0:.1f}".format(round(aligned_identity,1))) + '\t' + str("{0:.1f}".format(round(read_ins,1))) + '\t' + str("{0:.1f}".format(round(read_deln,1))) + '\t' + str("{0:.1f}".format(round(indels,1))) + '\t' + \
			#	str("{0:.0f}".format(round(perc_aligned,0))) + '\t' + alignment[1] + '\t' + str("{0:.0f}".format(round(alignment_string_length,0))) + '\t' + \
			#	chaining_score1 + '\t' + chaining_score2 + '\t' + aln_score + '\t' + str(ref_edit_distance) + '\t' + strand_orientation.strip('\n') + '\t' + aln_type + '\t' + str(x)  + '\t' +sam_flag + '\n'
			#	newfile.write(str(writ))
		except:
			if alignment[0].startswith('@'):
				#newfile3.write(line)
				badlines += 1
			elif not alignment[0].startswith('@'): #and alignment[1] != 4 and alignment[2] != '*':
			#	writ2 = alignment[0] + '\t'+ alignment[1] + '\t'+ alignment[2] + '\n'
			#	newfile2.write(str(writ2))
				noalign += 1
			else:   
				badlines += 1
		p+=1
		bar.update(p)
	bar.finish()
	k=0
	try:
		print('There were ' + str(badlines) + ' reads without alignments in the file.' + '\n' + 'Percentage of reads with alignments = ' + str("{0:.0f}".format(round(rdsaln/(noalign + rdsaln) * 100,0))))
	except:
		k+=1
	try:
		print('There were ' + str(badlines) + ' lines starting with @__ and ' + str(noalign) + ' lines with no alignments and are recorded.' + '\n' + 'Percentage of lines showing alignments = ' + str("{0:.0f}".format(round(rdsaln/(badlines + rdsaln + noalign) * 100,0))))
	except:
		k+=1	
	print('Print failure = ' + str(k))
	print('Finised, bye')
	newfile.close()
	newfile2.close()
	#newfile3.close()	

dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/capitata/rna_seq/Cc_E_5H/combined/minimap2_moy'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/eleftherios/minimap2'

#sam_analyse_bwa_mem(dr+'/'+'1.sam')
#sam_analyse_bwa_mem(dr+'/'+'C004_FailPass_2DnTemplate.contigs.sam')
#sam_analyse_bwa_mem(dr+'/'+'pilon.corrected.1.sam')
#sam_analyse_bwa_mem(dr+'/'+'pilon.corrected.2.sam')
#sam_analyse_bwa_mem(dr+'/'+'olivefly_supernova_v1p1_corrvalid_group_gt10_500000bc_74X_seed0.pseudohap.1.sam')
#sam_analyse_bwa_mem('/mnt/parallel_scratch_mp2_wipe_on_december_2017/ioannisr/banthony/analysis/nanopore/olivefly/assembly_annotation/illumina_alignment/olivefly_supernova_500000bc_74X/subsample_alignment/2.sam')
#sam_analyse_bwa_mem('/mnt/parallel_scratch_mp2_wipe_on_december_2017/ioannisr/banthony/analysis/nanopore/olivefly/combined/canu/polished/illumina_polished/quality_check/sam_analyse/pilon.corrected.2.sam')
#sam_analyse_bwa_mem('/mnt/parallel_scratch_mp2_wipe_on_december_2017/ioannisr/banthony/analysis/nanopore/olivefly/assembly_annotation/illumina_alignment/olivefly_supernova_500000bc_74X/sam_analyse/olivefly_supernova_v1p1_corrvalid_group_gt10_500000bc_74X_seed0.pseudohap.1.sam')
#sam_analyse_bwa_mem2('/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/bwa/tig00006001_pilon_pilon_tony-Bo.E.5H_all_rnaseq_combined_pass.trimmed_stranded1.sam')
#sam_analyse_graphmap_default(dr+'/'+'raw_reads_ercc_transcriptome_graphmap.sam')
#sam_analyse_graphmap_default(dr+'/'+'corrected_reads_ercc_transcriptome_graphmap.sam')
#maf_analyse_lastal('/gs/project/wst-164-ab/anthony/Nanopore_data/ecoli/ecoliTest/C002_05_6_50kb/nanonet/lastal/C002_05_6_50kb_nanonet_lastal.maf')
#sam_analyse_gmap_default('/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/capitata/rna_seq/gmap_all/GCF_000347755.2_Ccap_1.1_genomic_Cc.E.5H_all_rnaseq_combined_pass.trimmed_stranded.gmap.sam')
#sam_analyse_minimap2(dr+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs_canu_uncorrected_reads.names_2_minimap2.sam') #,dr+'/'+'chimera_reads.txt')
#sam_analyse_minimap2(dr+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs_canu_uncorrected_reads.names_2_lordec_corrected_minimap2.sam')
sam_analyse_minimap2(dr+'/'+'NW_013584145.1.85_contigs.sam')
#sam_analyse_minimap2(dr+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs_tapis_corrected_reads_for_isoforms_minimap.sam')
#sam_analyse_minimap2(dr+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs_sqanti_corrected_reads_for_isoforms_minimap.sam')

#sam_analyse_minimap2(dr+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs_Bo.E.5H_all_rnaseq_combined_pass.single_adapter.correctedReads_minimap2.sam')


#purse_cigar(dr+'/'+'E_coli_K12_NC000913.1_denovocomb.contigs2_graphmap.sam')

#sam_analyse_bwa_mem(dr+'/'+'C004_FailPass_2DnTemplate.contigs.sam')
#sam_analyse_bwa_mem(dr+'/'+'C004_FailPass_2DnTemplate.contigs_pilon1.sam')
#sam_analyse_bwa_mem(dr+'/'+'C004_FailPass_2DnTemplate.contigs_pilon2.sam')
#sam_analyse_bwa_mem(dr+'/'+'olivefly_supernova_v1p1.sam')