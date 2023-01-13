#!/usr/bin/env python3

import re, os, progressbar

def convertFasta(fasta):
	ofasta = open(fasta)
	rfasta = ofasta.read()
	seqPat = re.compile(r'([A-Z])\n([A-Z])')
	newFasta = seqPat.sub(r'\1\2', rfasta)	
	return(newFasta)
	
def break_seqs(fasta,fileformat):
	ofasta = open(fasta)
	rfasta = ofasta.readlines()	
	newfile = open(os.path.join(os.path.dirname(fasta), '_ed.'.join([os.path.basename(fasta).replace('.'+fileformat,''),fileformat])),'w')
	
	#if not os.path.exists(newfile):
	#	newfile = open(newfile, 'w')
	#else:
	#	print('File already exists\n')
	brk = 100000
	n=0
	x=0
	if rfasta[0].startswith('>') and not rfasta[2].startswith('>'):
		rfasta = convertFasta(fasta)
		rfasta = rfasta.strip().split('\n')
		bar = progressbar.ProgressBar(maxval=len(rfasta), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		bar.start()
		for line in rfasta:
			if not line.startswith('>'):
				seq = rfasta[rfasta.index(line)-1]
				length = len(line.strip())
				step = length/brk
				if '.' in str(step):
					k = int(str(step).split('.')[0])
				for i in range(0,k+1):
					newfile.write(seq.strip()+'_'+str(i)+'\n')
					if len(line[n:]) < brk:
						newfile.write(line[n:]+'\n')
					else:
						newfile.write(line[n:n+brk]+'\n')
					n+=brk
			n=0
			x+=1
			bar.update(x)
		bar.finish()
		newfile.close()
	else:
		bar = progressbar.ProgressBar(maxval=len(rfasta), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		bar.start()
		for line in rfasta:
			if not line.startswith('>'):
				seq = rfasta[rfasta.index(line)-1]
				length = len(line.strip())
				step = length/brk
				if '.' in str(step):
					k = int(str(step).split('.')[0])
				for i in range(0,k+1):
					newfile.write(seq.strip()+'_'+str(i)+'\n')
					if len(line[n:]) < brk:
						newfile.write(line[n:])
					else:
						newfile.write(line[n:n+brk]+'\n')
					n+=brk
			n=0
			x+=1
			bar.update(x)
		bar.finish()
		newfile.close()	


dr = '/media/abayega/Padlock_DT/analysis/nanopore/ecoli/C002-05-4/circleator/'
#break_seqs(dr+'test.fasta','fasta')
break_seqs(dr+'denovocomb.contigs.fasta','fasta')