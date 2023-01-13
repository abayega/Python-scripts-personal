#!/usr/bin/env python

def break_scaffolds_2_contigs(fasta, file_format):
	newfasta = open(os.path.dirname(fasta)+'/'+os.path.basename(fasta)+'_contigs', 'w')
	o = 0
	k = 0
	lenfasta = int(os.popen('wc -l %s' %(fasta)).read().split()[0])
	bar = progressbar.ProgressBar(maxval=lenfasta, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	seqPat = re.compile(r'(N+)')
	newFasta = seqPat.sub(r'\1\2', rfasta)
	while x<lenfasta:
		a=file1.readline()
		b=file1.readline().strip()
		if 'N' in b:
			Ns = b.count("N")


