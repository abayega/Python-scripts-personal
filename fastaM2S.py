#!/usr/bin/env python3

import re, os, progressbar

def convertFasta(fasta,fileformat):
	ofasta = open(fasta)
	rfasta = ofasta.read()
	#s.upper()
	print('Working')
	seqPat = re.compile(r'([A-Z])\n([A-Z])')
	newFasta = seqPat.sub(r'\1\2', rfasta)

	newfile_name = os.path.join(os.path.dirname(fasta), '_sl.'.join([os.path.basename(fasta).replace('.'+fileformat,''),fileformat]))
	if not os.path.exists(newfile_name):		
		newfile = open(newfile_name, 'w')
		newfile.write(newFasta)
		newfile.close()	
		print('Done')
		
def convertFasta2(fasta, fmt):
	#newfile_name = os.path.join(os.path.dirname(fasta), '_sl.'.join([os.path.basename(fasta).replace('.'+fileformat,''),fileformat]))
	newfile_name = os.path.dirname(fasta) + '/' + os.path.basename(fasta).replace('.'+fmt,'_sl.fasta')
	newfile = open(newfile_name, 'w')
	#if not os.path.exists(newfile_name):
		#newfile = open(newfile_name, 'w')
	reads = int(os.popen('grep ">" %s | wc -l' %(fasta)).read().split()[0])
	print('\nDetected ' + str(reads) + ' reads\n')
	length = int(os.popen('wc -l %s' %(fasta)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines\n')
	file1=open(fasta,'r')
	k=0
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while k < length:
		a=file1.readline()
		if not a.startswith('>'):
			newfile.write(a.strip())
		elif k==0 and a.startswith('>'):
			newfile.write(a)
		else:
			newfile.write('\n'+a)
		k+=1
		bar.update(k)
	bar.finish()
	newfile.close()

dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/ncbi_genome'
dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq/ERCC92'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/eleftherios/importin/dmel_importins'
dr = '/home/banthony/projects/rrg-ioannisr/banthony/analysis/covid19/genome'

convertFasta2(dr + '/' + 'MN908947.3.fasta', 'fasta')
#convertFasta2(dr + '/' + 'Bo_E_all_pass_edited.correctedReads.trimmed_stranded_choped.fasta', 'fasta')
#convertFasta2(dr + '/' + 'Bo_E_3H_C010_08_pass_edited.strand_no_strand_choped.fasta', 'fasta')
#convertFasta2(dr + '/' + 'Bo_E_5H_pass_edited.strand_no_strand_choped.fasta', 'fasta')
#convertFasta(dr+'/'+'Bo_E_5H_pass_edited_b4_canu_correct.choped.fasta','fasta')
#convertFasta(dr+'/'+'Bo_E_5H_pass_edited_canu_corrected.choped.fasta','fasta')
#convertFasta(dr+'/'+'Bo_E_5H_pass_edited_combined_corr_uncorr.choped.fasta','fasta')
#convertFasta(dr+'/'+'Bo_E_5H_pass_edited_uncorrected.choped.fasta','fasta')
