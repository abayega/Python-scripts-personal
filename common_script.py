#!/usr/bin/env python3

import re, os, sys, check_fasta, progressbar, shutil

def check_path_exists(path):
	if os.path.exists(path):
		sys.exit('\n...Path to write results exists, please remove and try again...\n')
		
def fasta_to_dict(fasta):
	#This function takes a fastq file with single line sequence and returns a dictionary with sequence name and corresponding read
	fastaNameSeq = {}	
	#openfasta = open(fasta, 'r')
	
	try:
		fasta_ed = check_fasta.check_fasta_fmt(fasta)
	except:
		if not os.path.exists(fasta_ed):
			sys.exit(-1)
			
	x = 0
	
	length = int(os.popen('wc -l %s' %(fasta_ed)).read().split()[0]) #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
	
	file1=open(fasta_ed,'r')
	bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length: #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
		name=file1.readline()
		sequence=file1.readline()
		fastaNameSeq.setdefault(name.strip()+'\n', sequence.strip()+'\n')	
		x+=2
		bar.update(x)
	bar.finish()
	return(fastaNameSeq)

def fastq_to_dict(fastq):
	#This function takes a fastq file with single line sequence and returns a dictionary with sequence name and corresponding read
	fastqNameSeq = {}	
	openfastq = open(fastq, 'r')
	
	x = 0
	
	length = int(os.popen('wc -l %s' %(fastq)).read().split()[0])
	file1=open(fastq,'r')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length:
		a=file1.readline()
		new_a=a.split()[0].strip('@')
		b=file1.readline()
		c=file1.readline()
		d=file1.readline()
		
		fastqNameSeq.setdefault(new_a, a+b+c+d.strip()+'\n')		
		x+=4
		bar.update(x)
	bar.finish()
	
	#print(fastqNameSeq)
	return(fastqNameSeq)

def get_reverse_complement(fasta):	
	complement = {'a':'T','c':'G','g':'C','t':'A','n':'N','A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
	seq=''
	for item in fasta[::-1]:
		seq=seq+complement[item]
	return seq
	
def get_metrichor_minus_nanonet(filetoread):
	fi = open(filetoread)
	lfi = fi.readlines()
	n = 0
	ufile = open(os.path.join(os.path.dirname(filetoread), '_'.join([os.path.basename(filetoread).rstrip('.xxx'),'met-nan.txt'])), 'w')
	
	for i in range(len(lfi)):
		if i%2 == 0:
			try:
				writ = int(lfi[i]) -  int(lfi[i+1])
				writ = str(writ) + '\n'
			except:
				n += 1
			ufile.write(writ)
	ufile.close()
	print(n)
	
def files_with_pattern():
	afiles = []
	files = os.listdir('/home/abayega/Documents/anthony/Nanopore_data/ecoliTest/ecoli50kb/downloads/fail')
	for i in files:
		if 'changed' in i:
			afiles = afiles + [i]
	print(len(afiles))
	
def remove_patern_from_file(fasta,pattern,fileformat):
	ofasta = open(fasta)
	rfasta = ofasta.read()
	seqPat = re.compile(r'%s' %(pattern))
	newFasta = seqPat.sub(r'', rfasta)
	
	newfile_name = os.path.join(os.path.dirname(fasta), '_1.'.join([os.path.basename(fasta).replace('.'+fileformat,''),fileformat]))
	if not os.path.exists(newfile_name):
		newfile = open(newfile_name, 'w')
		newfile.write(newFasta)
		newfile.close()
		
def find_pattern(gff3):
	newfile = open(gff3+'_genes', 'w')
	seqPat = re.compile(r'gene(\d+)')
	
	length = int(os.popen('wc -l %s' %(gff3)).read().split()[0])
	x = 0
	bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	
	for line in open(gff3):
		gene_ID = seqPat.search(line)
		if gene_ID:
			gene_ID = 'gene'+gene_ID.group(1)
			newfile.write(gene_ID+'\n')
		x += 1
		bar.update(x)
	bar.finish()
	newfile.close()
		
def find_gtf_geneID(gtf):
	newfile = open(gtf+'_genelist', 'w')
	ogtf = open(gtf)
	rgtf = ogtf.readlines()
	k = 0
	bar = progressbar.ProgressBar(maxval=len(rgtf), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in rgtf:
		linex = line.split()
		if 'gene_id' in linex:
			newfile.write(linex[linex.index('gene_id')+1].strip(';"') + '\n')
		k += 1
		bar.update(k)
	bar.finish()		
	newfile.close()
		
def cigar_to_length(file):
	newfile = open(file+'_exons', 'w')
	for cigar in open(file):
		#cigar = '23M1I15M3I1M1I8M1D3M1D9M1I17M1D7M1I16M1D18M1I21M1D9M2D10M'
		cigar_len = sum(int(x) for x in os.popen("echo %s | sed 's/[A-Z]/+/g'" %(cigar)).read().split()[0].strip('+').split('+'))
		newfile.write(cigar_len + '\n')
	newfile.close()
	
def get_ercc_mapping(ercc_names,rd_len_gc,er_len_gc,rd_er):
	#This script takes 4 files: a file with ercc names,a file with read name & the read length & read GC, a file with ercc read name & ercc read length & GC, and finally a file with
	#read name and the ercc read name it aligns to. The script goes over all names in file1 and returns for each the ercc read name, the ratio of read len:ercc len (that is observed:expected)
	#and the ratio of read GC:ercc GC.
	dirc = os.path.dirname(ercc_names)
	rercc_names = open(ercc_names)
	lercc_names = rercc_names.readlines()
	rrd_len_gc = open(rd_len_gc)
	lrd_len_gc = rrd_len_gc.readlines()
	rer_len_gc = open(er_len_gc)
	ler_len_gc = rer_len_gc.readlines()
	rrd_er = open(rd_er)
	lrd_er = rrd_er.readlines()
	newfile = open(dirc + '/' + 'ERCC_ratios', 'w')
	newfile.write('ercc_name\t'+'len_ratio\t'+'GC_ratio\t'+'read_len\t'+'read_GC\t'+'ercc_length\t'+'ercc_GC\t'+'\n')
	m = 0
	bar = progressbar.ProgressBar(maxval=len(lercc_names), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in lercc_names:
		xline = ([x for x in ler_len_gc if line.strip() in x][0])
		#print([xline])
		len1 = xline.split('\t')[1].strip()
		GC1 = xline.split('\t')[2].strip()
		ler_len_gc.remove(xline)
		os.system('grep %s %s | cut -f 1 > %s/line_mappings' %(line.strip(), rd_er, dirc))
		for eline in open(dirc+'/'+'line_mappings'):
			try:
				yline = ([x for x in lrd_len_gc if eline.strip() in x][0])
				len2 = yline.split('\t')[1].strip()
				GC2 = yline.split('\t')[2].strip()
				len_ratio = int(len2)/int(len1)
				GC_ratio = int(GC2)/int(GC1)
				lrd_len_gc.remove(yline)
				newfile.write(line.strip() + '\t' + str(len_ratio) + '\t' + str(GC_ratio) + '\t' + str(len2) + '\t' + str(GC2) + '\t' + str(len1) + '\t' + str(GC1) + '\n')
			except:
				print([eline])
		m += 1		
		bar.update(m)
	bar.finish()
	newfile.close()
	
def smthing(ercc_mappings,er_len_gc):
	dirc = os.path.dirname(ercc_mappings)
	newfile = open(dirc + '/' + 'ERCC_ratios_full', 'w')
	newfile.write('ercc_name\t'+'len_ratio\t'+'GC_ratio\t'+'read_len\t'+'read_GC\t'+'ercc_length\t'+'ercc_GC\t'+'\n')	
	rercc_mappings = open(ercc_mappings)
	lercc_mappings = rercc_mappings.readlines()
	rer_len_gc = open(er_len_gc)
	ler_len_gc = rer_len_gc.readlines()
	m = 0
	bar = progressbar.ProgressBar(maxval=len(lercc_mappings), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for i in ler_len_gc:
		ercc = i.split('\t')
		ercc_name = ercc[0]
		list1 = ([x for x in lercc_mappings if ercc_name.strip() in x])
		for e in list1:
			newfile.write(e.strip()+'\t'+ercc[1]+'\t'+ercc[2]) #+'\n')
			lercc_mappings.remove(e)
		m += 1		
		bar.update(m)
	bar.finish()
	newfile.close()

def dist_reads(reads_fasta):
	dirc = os.path.dirname(reads_fasta)
	base = os.path.basename(reads_fasta)
	nfile1 = open(dirc+'/1H/'+base.replace('.fasta','.1H.fasta'), 'w')
	nfile2 = open(dirc+'/2H/'+base.replace('.fasta','.2H.fasta'), 'w')
	nfile3 = open(dirc+'/3H/'+base.replace('.fasta','.3H.fasta'), 'w')
	nfile4 = open(dirc+'/4H/'+base.replace('.fasta','.4H.fasta'), 'w')
	nfile5 = open(dirc+'/5H/'+base.replace('.fasta','.5H.fasta'), 'w')
	nfile6 = open(dirc+'/6H/'+base.replace('.fasta','.6H.fasta'), 'w')
	nfile7 = open(dirc+'/female/'+base.replace('.fasta','.female.fasta'), 'w')
	nfile8 = open(dirc+'/male/'+base.replace('.fasta','.male.fasta'), 'w')
	nfile9 = open(dirc+'/'+'unknown_source_uncorrected.fasta', 'w')
	
	x=c=d=e=f=g=h=i=j=k = 0
	
	length = int(os.popen('wc -l %s' %(reads_fasta)).read().split()[0])
	print('\nWe are looping over ' + str(length/2) + ' reads\n')
	file1=open(reads_fasta,'r')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length:
		a=file1.readline()
		b=file1.readline()
		if "Bo_E_1" in a: #.split()[0]:
			nfile1.write(a+b)
			c += 1
		elif "Bo_E_2" in a: #.split()[0]:
			nfile2.write(a+b)
			d += 1
		elif "Bo_E_3" in a: #.split()[0]:
			nfile3.write(a+b)
			e += 1
		elif "Bo_E_4" in a: #.split()[0]:
			nfile4.write(a+b)
			f += 1
		elif "Bo_E_5" in a: #.split()[0]:
			nfile5.write(a+b)
			g += 1
		elif "Bo_E_6" in a: #.split()[0]:
			nfile6.write(a+b)
			h += 1
		elif "female" in a: #.split()[0]:
			nfile7.write(a+b)
			j += 1
		elif "male" in a: #.split()[0]:
			nfile8.write(a+b)
			k += 1
		else:
			nfile7.write(a+b)
			i += 1
			
		x+=2
		bar.update(x)
	bar.finish()
	if length/2 != sum([c,d,e,f,g,h,i,j,k]):
		print('The total lines in file does not equate number of reads distributed\n' + str(length/2) + ' ' + str(sum([c,d,e,f,g,h,i])))
	else:
		print('\nThe total lines in file is ' + str(length/2) + ' and total reads distributed is ' + str(sum([c,d,e,f,g,h,i,j,k])))
	print('There were %s, %s, %s, %s, %s, %s, %s, %s, %s corrected reads' %(str(c),str(d),str(e),str(f),str(g),str(h),str(j),str(k),str(i),))

def dist_reads2(sam_file):
	dirc = os.path.dirname(sam_file)
	base = os.path.basename(sam_file)
	nfile1 = open(dirc+'/'+base.replace('Bo.E.Heads','Bo_1H'), 'w') #Bo_E_1H_C010_10
	nfile2 = open(dirc+'/'+base.replace('Bo.E.Heads','Bo_2H'), 'w') #Bo_E_2H_C010_09
	nfile3 = open(dirc+'/'+base.replace('Bo.E.Heads','Bo_3H'), 'w') #Bo_E_3H_C010_08
	nfile4 = open(dirc+'/'+base.replace('Bo.E.Heads','Bo_4H'), 'w') #Bo_E_4H_C010_06
	nfile5 = open(dirc+'/'+base.replace('Bo.E.Heads','Bo_5H'), 'w')
	nfile6 = open(dirc+'/'+base.replace('Bo.E.Heads','Bo_6H'), 'w') #Bo_E_6H_C010_07
	nfile7 = open(dirc+'/'+base.replace('Bo.E.Heads','Heads_male'), 'w')
	nfile8 = open(dirc+'/'+base.replace('Bo.E.Heads','Heads_female'), 'w')
	nfile9 = open(dirc+'/'+'unknown_source_uncorrected.fasta', 'w')
	
	x=c=d=e=f=g=h=i=j=k=l = 0
	
	length = int(os.popen('wc -l %s' %(sam_file)).read().split()[0])
	print('\nWe are looping over ' + str(length) + ' lines\n')
	file1=open(sam_file,'r')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length:
		a=file1.readline()
		if "Bo_E_1" in a.split()[0]:
			nfile1.write(a)
			c += 1
		elif "Bo_E_2" in a.split()[0]:
			nfile2.write(a)
			d += 1
		elif "Bo_E_3" in a.split()[0]:
			nfile3.write(a)
			e += 1
		elif "Bo_E_4" in a.split()[0]:
			nfile4.write(a)
			f += 1
		elif "Bo_E_5" in a.split()[0]:
			nfile5.write(a)
			g += 1
		elif "Bo_E_6" in a.split()[0]:
			nfile6.write(a)
			h += 1
		elif "Bo_Heads_male" in a.split()[0]:
			nfile7.write(a)
			i += 1
		elif "Bo_Heads_female" in a.split()[0]:
			nfile8.write(a)
			j += 1
		elif a.startswith("@"):
			nfile1.write(a)
			nfile2.write(a)
			nfile3.write(a)
			nfile4.write(a)
			nfile5.write(a)
			nfile6.write(a)
			nfile7.write(a)
			nfile8.write(a)
			k += 1
		else:
			nfile9.write(a)
			l += 1
			
		x+=1
		bar.update(x)
	bar.finish()
	if length != sum([c,d,e,f,g,h,i,j,k,l]):
		print('\nThe total lines in file does not equate number of lines distributed...' + str(length) + ' vs ' + str(sum([c,d,e,f,g,h,i,j,k,l])))
	else:
		print('\nThe total lines in file is ' + str(length) + ' and total reads distributed is ' + str(sum([c,d,e,f,g,h,i,j,k,l])))
	print('\nThere were %s, %s, %s, %s, %s, %s, %s, %s, %s, %s lines, respectively...\n' %(str(c),str(d),str(e),str(f),str(g),str(h),str(i),str(j),str(k),str(l)))

def dist_reads3(sam_file):
	dirc = os.path.dirname(sam_file)
	base = os.path.basename(sam_file)
	nfile1 = open(dirc+'/'+base.replace('sam','1.sam'), 'w')
	nfile2 = open(dirc+'/'+base.replace('sam','2.sam'), 'w')
	nfile3 = open(dirc+'/'+base.replace('sam','3.sam'), 'w')
	nfile9 = open(dirc+'/'+'unknown_source_uncorrected.fasta', 'w')
	
	x=c=d=e=f=g=h=i=j=k=l = 0
	
	length = int(os.popen('wc -l %s' %(sam_file)).read().split()[0])
	print('\nWe are looping over ' + str(length) + ' lines\n')
	file1=open(sam_file,'r')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length:
		a=file1.readline()
		if a.startswith("@"):
			nfile1.write(a)
			nfile2.write(a)
			nfile3.write(a)
			k += 1
		elif not a.startswith("@") and x < 4281305:
			nfile1.write(a)
			d += 1
		elif not a.startswith("@") and x >= 4281305 and x < (4281305*2):
			nfile2.write(a)
			e += 1
		elif not a.startswith("@") and x >= (4281305*2):
			nfile3.write(a)
			f += 1
		else:
			nfile9.write(a)
			l += 1
			
		x+=1
		bar.update(x)
	bar.finish()
	if length != sum([c,d,e,f,g,h,i,j,k,l]):
		print('\nThe total lines in file does not equate number of lines distributed...' + str(length) + ' vs ' + str(sum([c,d,e,f,g,h,i,j,k,l])))
	else:
		print('\nThe total lines in file is ' + str(length) + ' and total reads distributed is ' + str(sum([c,d,e,f,g,h,i,j,k,l])))
	print('\nThere were %s, %s, %s, %s, %s, %s, %s, %s, %s, %s lines, respectively...\n' %(str(c),str(d),str(e),str(f),str(g),str(h),str(i),str(j),str(k),str(l)))
	
def dist_reads4(fasta):
	dirc = os.path.dirname(fasta)
	lopenfasta = fasta_to_dict(fasta)
	new_dic = {}
	bar_Pat = re.compile(r'C1H_(\w+)_FL')
	#for name,seq in lopenfasta.items():
	#	code = bar_Pat.search(name).group(1)
	#	cmd = 'echo "%s" >> "%s" ' %(name+seq, dirc+'/'+ code +'.fasta')
	#	os.system(cmd)
	for name,seq in lopenfasta.items():
		code = bar_Pat.search(name).group(1)
		if code not in new_dic.keys():
			new_dic.setdefault(code,[name+seq])
		else:
			new_dic[code].append(name+seq)
		#del lopenfasta[name]
	del lopenfasta
	for code_x,val in new_dic.items():
		newfile = open(dirc+'/'+code_x+'.fasta', 'w')
		for item_y in val:
			newfile.write(item_y)
		newfile.close()		
		
def add_rdlen_readname(reads_fasta):
	try:
		reads_fasta = check_fasta.check_fasta_fmt(reads_fasta)
	except:
		if not os.path.exists(reads_fasta):
			sys.exit(-1)
	nfile1 = open(reads_fasta+'_2', 'w')
	file1=open(reads_fasta,'r')
	x = 0
	length = int(os.popen('wc -l %s' %(reads_fasta)).read().split()[0])
	while x < length:
		a=file1.readline()
		b=file1.readline()
		#nfile1.write(a.split()[0]+str(len(b.strip())) + '\n' + b)
		nfile1.write(a.split()[0]+'|'+a.split()[1]+'|'+str(len(b.strip())) + '\n' + b)
		x += 2
	nfile1.close()
	
def add_rdlen(reads_txt):
	nfile1 = open(reads_txt+'_2', 'w')
	file1=open(reads_txt,'r')
	x=e=f=z = 0
	length = int(os.popen('wc -l %s' %(reads_txt)).read().split()[0])
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length:
		a=file1.readline()
		b = a.strip().split()		
		if x == 0:
			#nfile1.write('\t'.join([b[0],b[1],b[2],b[3],b[4],b[5],b[9],b[11],b[16],b[17]]) + '\t' + 'ref_len' + '\n') #minimap2
			nfile1.write('\t'.join([b[0],b[1],b[2],b[4],b[5],b[10],b[13],b[15],b[23],b[26]]) + '\t' + 'ref_len' + '\n') #gmap
			e += 1
		#elif b[17] != 'S':
		#	nfile1.write('\t'.join([b[0],b[1],b[2],b[3],b[4],b[5],b[9],b[11],b[16],b[17]]) + '\t' + b[1].split('|')[-1] + '\n')
		#	f += 1
		else:
			#nfile1.write('\t'.join([b[0],b[1],b[2],b[3],b[4],b[5],b[9],b[11],b[16],b[17]]) + '\t' + b[1].split('|')[-1] + '\n') #minimap2
			nfile1.write('\t'.join([b[0],b[1],b[2],b[3],b[4],b[9],b[12],b[14],b[22],b[25]]) + '\t' + b[1].split('|')[-1] + '\n') #gmap
			z += 1
		x += 1
		bar.update(x)
	bar.finish()
	print(str(length) + '\t' + str(e+f+z) + '\t' + str(z))
	nfile1.close()

def change_read_names(reads_fasta,gpd):
	#This script takes a fasta file and will manipulate the name of the file
	#The script also takes a gpd (genepred) file and will manipulate the second field of each line to change the contig name
	dirc = os.path.dirname(reads_fasta)
	nfile1 = open(reads_fasta+'_2', 'w')
	nfile2 = open(gpd+'_2', 'w')
	
	
	x=c=d=e=f=g=h=i = 0
	
	length = int(os.popen('wc -l %s' %(reads_fasta)).read().split()[0])
	print('We are looping over ' + str(length/2) + ' reads')
	file1=open(reads_fasta,'r')
	lfile1=file1.readlines()
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for i in lfile1:
		if i.startswith(">"):
			nfile1.write('_'.join(i.strip().split('.'))+'\n')
			c += 1
		else:
			nfile1.write(i)			
		x+=1
		bar.update(x)
	bar.finish()

	print('The sequences in file is ' + str(c))
	
	file2=open(gpd,'r')
	length = int(os.popen('wc -l %s' %(gpd)).read().split()[0])
	print('We are looping over ' + str(length) + ' gpd lines')
	x = 0
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length:
		a=file2.readline()
		a = a.strip().split()
		nfile2.write('\t'.join([a[0],'_'.join(a[1].split('.'))])+ '\t' + '\t'.join(a[2:])+'\n')
		d += 1			
		x+=1
		bar.update(x)
	bar.finish()
	print('The total lines in file is ' + str(d))
	
def manip_gff3(gff3):
	dirc = os.path.dirname(gff3)
	nfile1 = open(dirc+'/'+'b_oleae_genes_n_associated_proteins', 'w')
	#nfile1.write('chromosome\t' + 'transcript_name\t' + 'geneID\n')
	#nfile1.write('chromosome\t' + 'gene_name\t' + 'geneID\t' + 'gene_length\n') #If u go to NCBI Gene database u can search this gene symbol
	#nfile1.write('Chromosome\t' + 'Gene\t' + 'GeneID\t' + 'GeneBank\t' + 'Product\t' + 'NCBI_gene_symbol\n') #If u go to NCBI Gene database u can search this gene symbol
	
	c=d=e=f=g = 0
	items = {}
	length = int(os.popen('wc -l %s' %(gff3)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	file1=open(gff3,'r')
	lfile1=file1.readlines()
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for i in lfile1:
		ftrs = i.split('\t')
		'''
		if len(ftrs) > 3 and 'mRNA' in ftrs[2]:
			ftrs2 = i.split(";")
			product = [x for x in ftrs2 if x.startswith("product=")][0].split("=")[1]
			gene = [x for x in ftrs2 if x.startswith("Parent=")][0].split("=")[1]
			gene2 = [x for x in ftrs2 if x.startswith("gene=")][0].split("=")[1]
			Dbxref = [x for x in ftrs2 if x.startswith("Dbxref")]
			GeneIDGB = Dbxref[0].split(",")
			GeneID = GeneIDGB[0].split(":")[1]
			Genbank = GeneIDGB[1].split(":")[1]
			if not product.startswith("uncharacterized"):
				c += 1
				nfile1.write(ftrs[0]+'\t'+gene+'\t'+GeneID+'\t'+Genbank+'\t'+product+'\t'+gene2+'\n')
			elif product.startswith("uncharacterized"):
				e += 1
			else:
				f += 1
			g += 1
		'''
		#'''
		if len(ftrs) > 3 and 'gene' in ftrs[2]:
			GeneID = "."
			description = "."
			gene_length = int(ftrs[4]) - int(ftrs[3]) +1
			ftrs2 = ftrs[-1].split(";")
			gene_name=ftrs2[0].strip("ID=")
			Dbxref = [x for x in ftrs2 if x.startswith("Dbxref")]
			if len(Dbxref) > 0:
				GeneIDGB = Dbxref[0].split(",")
				GeneID = GeneIDGB[0].split(":")[1]
			pdescription = [x for x in ftrs2 if x.startswith("description")]
			if len(pdescription) > 0:
				description = pdescription[0].split("=")[1]
			if GeneID not in items.keys():
				#items.setdefault(GeneID,[gene_name+'\t'+GeneID])
				items.setdefault(GeneID,[gene_name,GeneID])
			#nfile1.write(gene_name+'\t'+GeneID+'\t'+description+'\n')
			#nfile1.write(ftrs[0]+'\t'+gene_name+'\t'+GeneID+'\t'+str(gene_length)+'\n')
			#nfile1.write(ftrs[0]+'\t'+ftrs[3]+'\t'+ftrs[4]+'\t'+gene_name+'\t'+'1000'+'\t'+'.'+'\t'+ftrs[3]+'\t'+ftrs[4]+'\t'+'240,59,32'+'\t'+'1'+	str(gene_length)+','+'\t'+'0,\n') #To make a bed file
			g += 1
		
		if len(ftrs) > 3 and 'CDS' in ftrs[2]:
			GeneID = "."
			Genbank = "."
			product = "."
			ftrs2 = ftrs[-1].split(";")
			Dbxref = [x for x in ftrs2 if x.startswith("Dbxref")]
			if len(Dbxref) > 0:
				GeneIDGB = Dbxref[0].split(",")
				for ix in GeneIDGB:
					if ix.startswith('Dbxref=GeneID'):
						GeneID = ix.split(":")[1]
					elif ix.startswith('Genbank:'):
						Genbank = ix.split(":")[1]
					elif ix.startswith('Dbxref=Genbank'):
						Genbank = ix.split(":")[1]
					elif ix.startswith('GeneID'):
						#GeneID = GeneIDGB[0].split(":")[1]
						GeneID = ix.split(":")[1]
			p_product = [x for x in ftrs2 if x.startswith("product")]
			if len(p_product) > 0:
				product = p_product[0].split("=")[1]
			if GeneID in items.keys():
				#if len(items[GeneID]) == 1:
				#	items[GeneID].append(Genbank+'\t'+product)
				if len(items[GeneID]) == 2:
					items[GeneID].append(product)
				if Genbank not in items[GeneID]:
					items[GeneID].append(Genbank)
				
			#nfile1.write(ftrs[0]+'\t'+gene_name+'\t'+GeneID+'\t'+str(gene_length)+'\n')
			#nfile1.write(ftrs[0]+'\t'+gene_name+'\t'+GeneID+'\t'+str(gene_length)+'\n')
			#nfile1.write(ftrs[0]+'\t'+ftrs[3]+'\t'+ftrs[4]+'\t'+gene_name+'\t'+'1000'+'\t'+'.'+'\t'+ftrs[3]+'\t'+ftrs[4]+'\t'+'240,59,32'+'\t'+'1'+	str(gene_length)+','+'\t'+'0,\n') #To make a bed file
			g += 1
		#''' 
		'''
		if len(ftrs) > 3 and 'exon' in ftrs[2]:			
			ftrs2 = i.split(";")
			#print(ftrs2)
			transcript = [x for x in ftrs2 if x.startswith("Parent=")][0].split("=")[1]
			Dbxref = [x for x in ftrs2 if x.startswith("Dbxref")]
			if len(Dbxref) > 0:
				GeneIDGB = Dbxref[0].split(",")
				GeneID = GeneIDGB[0].split(":")[1]
			else:
				GeneID = "."
			nfile1.write(ftrs[0]+'\t'+transcript+'\t'+GeneID+'\n')
			g += 1
		'''
		'''
		if len(ftrs) > 3 and 'tRNA' in ftrs[2]:
			gene_length = int(ftrs[4]) - int(ftrs[3]) +1
			ftrs2 = ftrs[-1].split(";")
			rna_name=ftrs2[0].strip("ID=")
			#print(ftrs2)
			Dbxref = [x for x in ftrs2 if x.startswith("Dbxref=")]
			if len(Dbxref) > 0:
				GeneID = Dbxref[0].strip("Dbxref=GeneID:")
			else:
				GeneID = "."
			gene_name = [x for x in ftrs2 if x.startswith("Parent=")]
			if len(gene_name) > 0:
				gene_name = gene_name[0].split("=")[1]
			else:
				gene_name = "."
			#nfile1.write(ftrs[0]+'\t'+gene_name+'\t'+GeneID+'\t'+str(gene_length)+'\n')
			nfile1.write(ftrs[0]+'\t'+ftrs[3]+'\t'+ftrs[4]+'\t'+gene_name+'_'+rna_name+'_'+GeneID+'\t'+'1000'+'\t'+'.'+'\t'+ftrs[3]+'\t'+ftrs[4]+'\t'+'240,59,32'+'\t'+'1'+	str(gene_length)+','+'\t'+'0,\n') #To make a bed file
			g += 1
		'''
		d+=1
		bar.update(d)
	bar.finish()
	for gene,itms in items.items():
		if len(itms)==2:
			itms += ['.\t.']
			#'\t'.join([itms[0],itms[1],','.joint(itms[2:])])
		nfile1.write('\t'.join([itms[0],itms[1],itms[2],','.join(itms[3:])])+'\n')
	nfile1.close()

	#print('Characterised features recorded for ' + str(c) + ' genes.')
	#print('Uncharacterised = ' + str(e))
	#print('Unknown =  ' + str(f))
	#print(g)
	
def manip_gff3_2(gff3):
	dirc = os.path.dirname(gff3)
	nfile1 = open(dirc+'/'+'b_oleae_genes_n_associated_proteins', 'w')
	nfile2 = open(dirc+'/'+'b_oleae_genes_n_associated_proteins', 'w')
	c=d=e=f=g = 0
	items = {}
	length = int(os.popen('wc -l %s' %(gff3)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	file1=open(gff3,'r')
	lfile1=file1.readlines()
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line_x in lfile1:
		ftrs = line_x.split('\t')
		
		if len(ftrs) > 3 and 'mRNA' in ftrs[2]:
			rnaIDPat = re.compile(r'ID=(\w+);')
			gene_IDPat = re.compile(r'Parent=(\w+);')
			Genbank_IDPat = re.compile(r'Genbank:(\w+.\w?);')
			trans_IDPat = re.compile(r'transcript_id=(\w+.\w?)')
			
			rnaID = gene_ID = Genbank_ID = trans_ID = "."
			try:
				rnaID = rnaIDPat.search(line).group(1)
				gene_ID = gene_IDPat.search(line).group(1)
				Genbank_ID = Genbank_IDPat.search(line).group(1)
				trans_ID = trans_IDPat.search(line).group(1)
			except:
				nfile2.write(line_x)
			if 	Genbank_ID not in items.keys():
				items.setdefault(Genbank_ID,':'.join([rnaID,gene_ID,Genbank_ID,trans_ID]))			
		d+=1
		bar.update(d)
	bar.finish()
	
	for gene,itms in items.items():
		if len(itms)==2:
			itms += ['.\t.']
			#'\t'.join([itms[0],itms[1],','.joint(itms[2:])])
		nfile1.write('\t'.join([itms[0],itms[1],itms[2],','.join(itms[3:])])+'\n')
	nfile1.close(); nfile2.close()
	
def gff3_to_gtf(gff3):
	#This script takes a gff3 and will try to turn it into a gtf usable by /home/banthony/software/gtex-pipeline-master/gene_model/collapse_annotation_original.py
	#which collapses all isoforms into a single transcript
	dirc = os.path.dirname(gff3)
	name = os.path.basename(gff3)
	nfile1 = open(dirc+'/'+name+'.gencode.gtf', 'w')
	#nfile2 = open(dirc+'/'+'gtf_new2', 'w')
	exon_num2 = 0
	a=c=d=e=f=g = 0
	items = {}
	length = int(os.popen('wc -l %s' %(gff3)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for i in open(gff3, 'r'):
		#nfile2.write(i)
		ftrs = i.strip().split('\t')
		if len(ftrs) > 3 and 'gene' in ftrs[2]:
			GeneID = "."
			description = "."
			gName='.'
			level = "NA"
			gene_length = int(ftrs[4]) - int(ftrs[3]) +1
			ftrs2 = ftrs[-1].split(";")
			gene_id=ftrs2[0].strip("ID=")
			Name=[x for x in ftrs2 if x.startswith("Name")]
			Dbxref = [x for x in ftrs2 if x.startswith("Dbxref")]
			try:
				gene_type = [x for x in ftrs2 if x.startswith("gene_biotype=")][0].split("=")[1].strip()
			except:
				gene_type = "unknown"
			if len(Name) > 0:
				gName=Name[0].split("=")[1]
			if len(Dbxref) > 0:
				GeneIDGB = Dbxref[0].split(",")
				GeneID = GeneIDGB[0].split(":")[1]
			pdescription = [x for x in ftrs2 if x.startswith("description")]
			if len(pdescription) > 0:
				description = pdescription[0].split("=")[1]
			if GeneID not in items.keys():
				items.setdefault(GeneID,[gene_id,gName,gene_type])
			g += 1
			nfile1.write('\t'.join(ftrs[:8])+'\t'+'gene_id "%s"; gene_name "%s"; gene_type "%s"; level "%s";\n' %(gene_id,gName,gene_type,level))
			
		if len(ftrs) > 3 and 'gene' in ftrs[2] and "pseudo=true" in i:
			transcript_id = "pseudo_"+str(a+1)
			a+=1
		if len(ftrs) > 3 and "pseudo=true" not in i:
			if 'mRNA' in ftrs[2] or 'ncRNA' in ftrs[2] or 'transcript' in ftrs[2] or 'tRNA' in ftrs[2] or 'rRNA' in ftrs[2]:
				ftrs2 = ftrs[-1].split(";")
				transcript_id=ftrs2[0].strip("ID=")
				transcript_type = "FU"
				try:
					product = [x for x in ftrs2 if x.startswith("product=")][0].split("=")[1]
					Parent = [x for x in ftrs2 if x.startswith("Parent=")][0].split("=")[1]
					gene2 = [x for x in ftrs2 if x.startswith("gene=")][0].split("=")[1]
					Dbxref = [x for x in ftrs2 if x.startswith("Dbxref")]
					GeneIDGB = Dbxref[0].split(",")
					GeneID = GeneIDGB[0].split(":")[1]
					transcript_name = [x for x in ftrs2 if x.startswith("Name=")][0].split("=")[1]
				except:
					Parent = transcript_id
					GeneID = transcript_id
					items.setdefault(GeneID,[transcript_id,transcript_id,product])
					transcript_name = "unknown"
				if not product.startswith("uncharacterized"):
					c += 1
				elif product.startswith("uncharacterized"):
					e += 1
				else:
					f += 1
				g += 1
				parent = items[GeneID]
				gene_id = parent[0]
				gName = parent[1]
				gene_type = parent[2]
				nfile1.write('\t'.join(ftrs[:8]).replace(ftrs[2],'transcript')+'\t'+'gene_id "%s"; gene_name "%s"; transcript_id "%s"; transcript_name "%s"; transcript_type "%s"; gene_type "%s"; level "%s";\n' 
							 %(gene_id,gName,transcript_id,transcript_name,transcript_type,gene_type,level))
			
		if len(ftrs) > 3 and 'exon' in ftrs[2]:			
			ftrs2 = ftrs[-1].split(";")
			exon_id=ftrs2[0].strip("ID=")
			exon_num=exon_id.strip('id')
			if not exon_num.isnumeric():
				exon_num = int(exon_num2)+1
			transcript = [x for x in ftrs2 if x.startswith("Parent=")][0].split("=")[1]
			Dbxref = [x for x in ftrs2 if x.startswith("Dbxref")]
			if len(Dbxref) > 0:
				GeneIDGB = Dbxref[0].split(",")
				GeneID = GeneIDGB[0].split(":")[1]
			else:
				GeneID = transcript
			try:
				parent = items[GeneID]
				gene_id = parent[0]
				gName = parent[1]
				gene_type = parent[2]
				nfile1.write('\t'.join(ftrs[:8])+'\t'+'gene_id "%s"; gene_name "%s"; transcript_id "%s"; exon_id "%s"; exon_number %s; gene_type "%s"; level "%s";\n'
							 %(gene_id,gName,transcript_id,exon_id,str(exon_num),gene_type,level))
			except:
				g += 1
			exon_num2 = exon_num
		d+=1
		bar.update(d)
	bar.finish()
	nfile1.close()
	#nfile2.close()
	
def gff3_to_gtf_ftr(gff3):
	#This script takes a gff3, pick mRNA genes, and will try to turn that into a gtf usable by /home/banthony/software/gtex-pipeline-master/gene_model/collapse_annotation_original.py
	#which collapses all isoforms into a single transcript
	dirc = os.path.dirname(gff3)
	name = os.path.basename(gff3)
	nfile1 = open(dirc+'/'+name+'.gencode.gtf', 'w')
	#nfile2 = open(dirc+'/'+'gtf_new2', 'w')
	exon_num2 = 0
	a=c=d=e=f=g = 0
	items = {}
	mRNA_genes = []
	length = int(os.popen('wc -l %s' %(gff3)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	#Go through the gff3 to record all gene IDs for mRNA genes
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for i in open(gff3, 'r'):
		ftrs = i.strip().split('\t')
		if len(ftrs) > 3 and 'mRNA' in ftrs[2]:
			ftrs2 = ftrs[-1].split(";")
			Dbxref = [x for x in ftrs2 if x.startswith("Dbxref")]
			GeneIDGB = Dbxref[0].split(",")
			GeneID = GeneIDGB[0].split(":")[1]
			mRNA_genes+=[GeneID]
			d+=1
		bar.update(d)
	bar.finish()
	d=0
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for i in open(gff3, 'r'):
		ftrs = i.strip().split('\t')
		if len(ftrs) > 3 and 'gene' in ftrs[2] and "protein_coding" in ftrs[-1]:
			GeneID = "."
			description = "."
			gName='.'
			level = "NA"
			gene_length = int(ftrs[4]) - int(ftrs[3]) +1
			ftrs2 = ftrs[-1].split(";")
			gene_id=ftrs2[0].strip("ID=")
			Name=[x for x in ftrs2 if x.startswith("Name")]
			Dbxref = [x for x in ftrs2 if x.startswith("Dbxref")]
			try:
				gene_type = [x for x in ftrs2 if x.startswith("gene_biotype=")][0].split("=")[1].strip()
			except:
				gene_type = "unknown"
			if len(Name) > 0:
				gName=Name[0].split("=")[1]
			if len(Dbxref) > 0:
				GeneIDGB = Dbxref[0].split(",")
				GeneID = GeneIDGB[0].split(":")[1]
			pdescription = [x for x in ftrs2 if x.startswith("description")]
			if len(pdescription) > 0:
				description = pdescription[0].split("=")[1]
			if GeneID not in items.keys():
				items.setdefault(GeneID,[gene_id,gName,gene_type])
			g += 1
			nfile1.write('\t'.join(ftrs[:8])+'\t'+'gene_id "%s"; gene_name "%s"; gene_type "%s"; level "%s";\n' %(gene_id,gName,gene_type,level))
			
		if len(ftrs) > 3 and "pseudo=true" not in i:
			if 'mRNA' in ftrs[2]:
				ftrs2 = ftrs[-1].split(";")
				transcript_id=ftrs2[0].strip("ID=")
				transcript_type = "FU"
				try:
					product = [x for x in ftrs2 if x.startswith("product=")][0].split("=")[1]
					Parent = [x for x in ftrs2 if x.startswith("Parent=")][0].split("=")[1]
					gene2 = [x for x in ftrs2 if x.startswith("gene=")][0].split("=")[1]
					Dbxref = [x for x in ftrs2 if x.startswith("Dbxref")]
					GeneIDGB = Dbxref[0].split(",")
					GeneID = GeneIDGB[0].split(":")[1]
					transcript_name = [x for x in ftrs2 if x.startswith("Name=")][0].split("=")[1]
				except:
					Parent = transcript_id
					GeneID = transcript_id
					items.setdefault(GeneID,[transcript_id,transcript_id,product])
					transcript_name = "unknown"
				if not product.startswith("uncharacterized"):
					c += 1
				elif product.startswith("uncharacterized"):
					e += 1
				else:
					f += 1
				g += 1
				parent = items[GeneID]
				gene_id = parent[0]
				gName = parent[1]
				gene_type = parent[2]
				if GeneID in mRNA_genes:
					nfile1.write('\t'.join(ftrs[:8]).replace(ftrs[2],'transcript')+'\t'+'gene_id "%s"; gene_name "%s"; transcript_id "%s"; transcript_name "%s";transcript_type "%s"; gene_type "%s"; level "%s";\n' %(gene_id,gName,transcript_id,transcript_name,transcript_type,gene_type,level))
			
		if len(ftrs) > 3 and 'exon' in ftrs[2] :		
			ftrs2 = ftrs[-1].split(";")
			exon_id=ftrs2[0].strip("ID=")
			exon_num=exon_id.strip('id')
			if not exon_num.isnumeric():
				exon_num = int(exon_num2)+1
			transcript = [x for x in ftrs2 if x.startswith("Parent=")][0].split("=")[1]
			Dbxref = [x for x in ftrs2 if x.startswith("Dbxref")]
			if len(Dbxref) > 0:
				GeneIDGB = Dbxref[0].split(",")
				GeneID = GeneIDGB[0].split(":")[1]
			else:
				GeneID = transcript
			try:
				parent = items[GeneID]
				gene_id = parent[0]
				gName = parent[1]
				gene_type = parent[2]
				if GeneID in mRNA_genes:
					nfile1.write('\t'.join(ftrs[:8])+'\t'+'gene_id "%s"; gene_name "%s"; transcript_id "%s"; exon_id "%s"; exon_number %s; gene_type "%s"; level "%s";\n' %(gene_id,gName,transcript_id,exon_id,str(exon_num),gene_type,level))
			except:
				g += 1
			exon_num2 = exon_num
		d+=1
		bar.update(d)
	bar.finish()
	nfile1.close()

def find_uniprot(uniprot_headers, proteins):
	gff3 = uniprot_headers
	dirc = os.path.dirname(gff3)
	nfile1 = open(dirc+'/'+'uniprot_annnotations', 'w')
	nfile1.write('Entry\t' + 'Entry_name\t' + 'Protein_name\t' + 'Organism\n') #If u go to NCBI Gene database u can search this gene symbol
	
	c=d=e=f=g = 0
	
	length = int(os.popen('wc -l %s' %(proteins)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	file1=open(gff3,'r')
	file2=open(proteins)
	lfile1=file1.readlines()
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(proteins):
		for i in lfile1:
			if line.strip() in i:
				entry = i.split("|")[1]
				entry_name = i.split("|")[2].split()[0]
				organism = i.split("OS=")[1].split("OX=")[0].strip()
				nfile1.write(entry+'\t'+entry_name+'\t'+line.strip()+'\t'+organism+'\n')
				lfile1.remove(i)
				c += 1
			#g += 1
		d+=1		
		bar.update(d)
	bar.finish()
	nfile1.close()

	print('Found ' + str(c) + ' Proteins.')

def find_uniprot2(uniprot_headers, proteins):
	dirc = os.path.dirname(uniprot_headers)
	nfile1 = open(dirc+'/'+'uniprot_trembl_annnotations_exact_match_removing-like', 'w')
	#nfile1x = open(dirc+'/'+'uniprot_trembl_annnotations_grep_match_headers_formated', 'w')
	nfile2 = open(dirc+'/'+'no_uniprot_trembl_annnotations_exact_match_after_removing-like', 'w')
	#nfile3 = open(dirc+'/'+'bad_format_grep', 'w')
	nfile1.write('Entry\t' + 'Entry_name\t' + 'Protein_name\t' + 'Organism\n')
	#nfile1x.write('Entry\t' + 'Entry_name\t' + 'Protein_name\t' + 'Organism\n') #If u go to NCBI Gene database u can search this gene symbol
	#nfile1.write('Query\t' +'Entry\t' + 'Entry_name\t' + 'Protein_name\t' + 'Organism\n') #If u go to NCBI Gene database u can search this gene symbol
	
	#'''#This code can be uncommented to create the necessary dictionary
	c=d=e=f=g = 0
	dictn = {}
	
	length = int(os.popen('wc -l %s' %(uniprot_headers)).read().split()[0])
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	
	for i in open(uniprot_headers):
		entry = i.split("|")[1]
		entry_name = i.split("|")[2].split()[0]
		organism = i.split("OS=")[1].split("OX=")[0].strip()
		protein_name = i.split(entry_name)[1].split("OS=")[0].strip()
		protein_name2 = protein_name.lower().split(",")[0] #added .split(",")[0] in order to remove the commas that might cause issues
		dictn.setdefault(entry+'\t'+entry_name+'\t'+protein_name+'\t'+organism+'\n',protein_name2) #This is incredibly powerful to get all hits
		#dictn.setdefault(protein_name2.lower(), entry+'\t'+entry_name+'\t'+protein_name+'\t'+organism+'\n')
		#dictn.setdefault(entry, entry+'\t'+entry_name+'\t'+protein_name+'\t'+organism+'\n')
		d+=1		
		bar.update(d)
	bar.finish()
	print(str(d)+' proteins added to dictionary\n')
	#'''
	
	''' #This code probably only runs slower but yeilds no more information
	length = int(os.popen('wc -l %s' %(proteins)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	
	for line in open(proteins):
		a = line.strip().lower()
		if len(a) > 6:
			for x,y in dictn.items():
				#if a in y.split('\t')[2].strip().lower(): #To get all ~ matches
				if a == y.split('\t')[2].strip().lower(): #To get exact matches
					nfile1.write(line.strip()+'\t'+y)
					g += 1
					c += 1
		if g == 0:
			nfile2.write(line)
		else:
			f += 1
			g = 0
		e+=1		
		bar.update(e)
	bar.finish()
	nfile1.close()
	nfile2.close()
	'''
	
	''' This code is to perform a grep so that you get unspecific matches
	c=d=e=f=g = 0
	file = dirc+'/'+'uniprot_trembl_annnotations_grep_match_headers'
	length = int(os.popen('wc -l %s' %(proteins)).read().split()[0])
	print('Starting to grep ' + str(length) + ' proteins...\n')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(proteins):
		if len(line.strip()) > 6:
			try:
				os.system('grep -i "%s" %s >> %s' %(line.strip(),uniprot_headers,file))
			except:
				nfile3.write(line)
			g += 1
			bar.update(g)
	bar.finish()
	nfile3.close()
	print('greping done. Now starting to write proteins in good format...\n')
	#file1 = dirc+'/'+'uniprot_trembl_annnotations_exact_match_headers'
	#file1 = dirc+'/'+'uniprot_sprot_annnotations_exact_match_headers'
	length = int(os.popen('wc -l %s' %(file)).read().split()[0])
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for i in open(file):
		entry = i.split("|")[1]
		entry_name = i.split("|")[2].split()[0]
		organism = i.split("OS=")[1].split("OX=")[0].strip()
		protein_name = i.split(entry_name)[1].split("OS=")[0].strip()
		nfile1x.write(entry+'\t'+entry_name+'\t'+protein_name+'\t'+organism+'\n') 
		d+=1		
		bar.update(d)
	bar.finish()
	nfile1x.close()
	print('Formatted and saved ' + str(d) + ' proteins homologues')
	'''
	
	#for line in open(proteins):
	#	if line.strip().lower() in dictn.keys():
	#		nfile1.write(dictn[line.strip().lower()])
	#		c += 1
	#	else:
	#		nfile2.write(line)
	#	e+=1		
	#	bar.update(e)
	#bar.finish()
	#nfile1.close()
	#nfile2.close()
	
	''' #This code can be uncommented
	dictn_keys = list(dictn.keys())
	for line in open(proteins):
		k = [i for i, x in enumerate(list(dictn.values())) if x == line.strip().lower()]
		if len(k) > 0:
			for item in k:
				nfile1.write(dictn_keys[item])
				c += 1
				#list(mydict.keys())[list(mydict.values()).index(16)]
			f += 1
		else:
			nfile2.write(line)
			e += 1
		e+=1		
		bar.update(e)
	bar.finish()
	nfile1.close()
	nfile2.close()

	print('Found ' + str(f) + ' Proteins.')
	print('Missing ' + str(e) + ' Proteins.')
	print('Total isoforms are ' + str(c))
	'''
	
	#'''
	c=d=e=f=g=h=j=k=l=n=m = 0	
	length = int(os.popen('wc -l %s' %(proteins)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	
	for line in open(proteins):
		a = line.strip().lower()
		#if len(a) > 6:
		for x,y in dictn.items():
			#if a in y.split('\t')[2].strip().lower(): #To get all ~ matches
			if a == y: #To get exact matches
				nfile1.write(x)
				k += 1
				g += 1
				c += 1
			elif a.replace("-like",'') == y: # a.rstrip("-like") == y: #To get exact matches otherwise not matching
				nfile1.write(x)
				g += 1
				c += 1
				h += 1
			elif a + "-like protein" == y:
				g += 1
				n += 1
		if g == 0:
			nfile2.write(line)
			e += 1
		else:
			f += 1
			g = 0
		if h > 0:
			j += 1
			h = 0
		if k > 0:
			l += 1
			k = 0
		if n > 0:
			m += 1
			n = 0
			
		d+=1		
		bar.update(d)
	bar.finish()
	nfile1.close()
	nfile2.close()

	print('Found ' + str(f) + ' Proteins.')
	print('Missing ' + str(e) + ' Proteins.')
	print('Total isoforms are ' + str(c))
	print('Exact matches are ' + str(l) + ' Proteins.')
	print('Proteins matching due to -like are ' + str(j))
	print(m)
	#'''
	
def correct_files(fasta):
	dirc = os.path.dirname(fasta)
	basename = os.path.basename(fasta)
	newfile = open(dirc+'/'+basename.replace('.fasta','_new.fasta'),'w')
	
	file1=open(fasta,'r')
	
	a=b=c=d=e=f=g=h=j=k=l=n=m=x = 0
	
	length = int(os.popen('wc -l %s' %(fasta)).read().split()[0])
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x<length:
		name=file1.readline().strip()
		name2=file1.readline().strip()
		seq=file1.readline().strip()
		
		if name2.startswith('.'):
			newfile.write(name+name2+'\n'+seq+'\n')
		x += 3
		a += 1
		bar.update(a)
	bar.finish()
	newfile.close()
	
def correct_files2(fasta):
	dirc = os.path.dirname(fasta)
	basename = os.path.basename(fasta)
	newfile = open(dirc+'/'+basename.replace('.fasta','_new.fasta'),'w')
	
	file1=open(fasta,'r')
	
	a=b=c=d=e=f=g=h=j=k=l=n=m=x = 0
	
	length = int(os.popen('wc -l %s' %(fasta)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x<length:
		name=file1.readline().strip()
		
		if name.startswith('>'):
			k = name
			newfile.write(k.split(".")[0].split()[0]+"."+".".join(k.split(".")[1:])+'\n')
			b += 1
		else:
			newfile.write(name+'\n')
			c += 1
		x += 1
		a += 1
		bar.update(a)
	bar.finish()
	newfile.close()
	print('Found ' + str(b) + ' sequence names, and ' + str(c) + ' sequence lines')
	
def correct_files3(fasta):
	dirc = os.path.dirname(fasta)
	basename = os.path.basename(fasta)
	newfile = open(dirc+'/'+basename.replace('.fasta','_new.fasta'),'w')
	
	file1=open(fasta,'r')
	
	a=b=c=d=e=f=g=h=j=k=l=n=m=x = 0
	m = 100
	b = 1
	length = int(os.popen('wc -l %s' %(fasta)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x<length:
		name=file1.readline().strip()
		
		if name.startswith('>'):
			k = name
			newfile.write("PB."+str(m+b)+"."+str(b)+"|"+name.lstrip(">")+"|blah"+'\n')
			b += 1
		else:
			newfile.write(name+'\n')
			c += 1
		x += 1
		a += 1
		bar.update(a)
	bar.finish()
	newfile.close()
	print('Found ' + str(b) + ' sequence names, and ' + str(c) + ' sequence lines')
	
def remove_bad_samlines(picard, samfile):
	dirc = os.path.dirname(samfile)
	basename = os.path.basename(samfile)
	newfile = open(dirc+'/'+basename.replace('.sam','_new.sam'),'w')
	newfile2 = open(dirc+'/'+basename.replace('.sam','_badlines.sam'),'w')
	a=b=c=d=e=f=g=h=j=k=l=n=m=x = 0
	samdic = {}
	file1=open(samfile,'r')
	length = int(os.popen('wc -l %s' %(samfile)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x<length:
		name=file1.readline()
		samdic.setdefault(name, name.split()[0])		
		x+=1
		bar.update(x)
	bar.finish()
	print("Finished creating sam dictionary, moving to removing bad lines...\n")
		
	a=b=c=d=e=f=g=h=j=k=l=n=m=x = 0
	
	file1=open(picard,'r')
	lines = []
	length = int(os.popen('wc -l %s' %(picard)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x<length:
		name=file1.readline().strip()
		name = name.split()[5].stirp(',')
		lines = lines + [name]
		x += 1
		bar.update(x)
	bar.finish()
	print("Finished getting names...")
	
	bar = progressbar.ProgressBar(maxval=len(samdic.values()), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for x,y in samdic.items():
		if y.strip() not in lines:
			newfile.write(x)
		elif y.strip() in lines:
			newfile2.write(x)
		a += 1	
		bar.update(a)
	bar.finish()
	newfile.close()
	print('Found ' + str(b) + ' sequence names, and ' + str(c) + ' sequence lines')

def remove_bad_samlines2(picard, samfile):
	dirc = os.path.dirname(samfile)
	basename = os.path.basename(samfile)
	newfile = open(dirc+'/'+basename.replace('.sam','_new.sam'),'w')
	newfile2 = open(dirc+'/'+basename.replace('.sam','_badlines.sam'),'w')
	a=b=c=d=e=f=g=h=j=k=l=n=m=x = 0
	
	file1=open(picard,'r')
	lines = []
	length = int(os.popen('wc -l %s' %(picard)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while a<length:
		name=file1.readline().strip()
		name = name.split()
		if len(name) > 6:
			name = name[5].strip(',')
			lines = lines + [name]
		a += 1
		bar.update(a)
	bar.finish()
	print("Finished getting names...")
	
	file1=open(samfile,'r')
	length = int(os.popen('wc -l %s' %(samfile)).read().split()[0])
	print('We are looping over ' + str(length) + ' sam-file lines')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x<length:
		name=file1.readline()
		name1 = name.split()
		if len(name1) > 6 and name1[0] not in lines:
			newfile.write(name)
			b += 1
		elif len(name1) > 6 and name1[0] in lines:
			newfile2.write(name)
			c += 1
		else:
			newfile.write(name)
		x+=1
		bar.update(x)
	bar.finish()
	print('Found ' + str(b) + ' good sam names, and ' + str(c) + ' bad sam lines')
		
def check_fasta_fmt(infile):
	import fastaM2S
	#Check fast
	length1 = int(os.popen('head -10000 %s | grep ">" | wc -l' %(infile)).read().split()[0])
	length2 = int(os.popen('head -10000 %s | grep -v ">" | wc -l' %(infile)).read().split()[0])
	if length1 == length2:
		print('\nFasta format is good...\n')
	else:
		print('\nFasta format is bad, creating a single line format, delete file if not needed...\n')
		fastaM2S.convertFasta2(infile)
		infile = path + '/' + os.path.basename(infile).replace('.fasta','_sl.fasta')
	
def change_read_names_to_pacbio(fasta):
	from Bio import SeqIO
	print('\Working on ' + fasta)
	try:
		fasta = check_fasta.check_fasta_fmt(fasta)
	except:
		if not os.path.exists(fasta):
			sys.exit(-1)
	length = int(os.popen('wc -l %s' %(fasta)).read().split()[0])/2 + 1
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	f = open(fasta.replace('.fasta','.modified.fasta'), 'w')
	i = 1
	for r in SeqIO.parse(open(fasta),'fasta'):
		f.write(">i{0}_HQ_Nanopore|cb{1}/f2p0/{2} {3}\n".format(i,i, len(r.seq), r.id))
		f.write("{0}\n".format(r.seq))
		i += 1
		bar.update(i)
	f.close()
	bar.finish()
	
def make_cluster_report_csv(collapsed_group_txt):
	dirc = os.path.dirname(collapsed_group_txt)
	basename = os.path.basename(collapsed_group_txt)
	newfile = open(dirc+'/'+'cluster_report.csv','w')
	
	a=b=c=d=e=f=g=h=j=k=l=n=m=x = 0
	
	file2 = open(collapsed_group_txt, 'r')
	length = int(os.popen('wc -l %s' %(collapsed_group_txt)).read().split()[0])
	print('\nWe are looping over ' + str(length) + ' groups\n')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x<length:
		group = file2.readline().strip()
		split_group = group.split()
		group_name = split_group[0]
		reads = split_group[1].split(',')
		for read in reads:
			newfile.write(group_name + '|' + read.split('|')[1] + '_CCS,FL\n')
		x += 1
		bar.update(x)
	bar.finish()
	print('\nDone, bye \n')
	
def write_header(newfile, total_reads):
	f = newfile
	f.write("#\n")
	f.write("# -----------------\n")
	f.write("# Field explanation\n")
	f.write("# -----------------\n")
	f.write("# count_fl: Number of associated FL reads\n")
	#f.write("# count_nfl: Number of associated FL + unique nFL reads\n")
	#f.write("# count_nfl_amb: Number of associated FL + unique nFL + weighted ambiguous nFL reads\n")
	f.write("# norm_fl: count_fl / total number of FL reads\n")
	#f.write("# norm_nfl: count_nfl / total number of FL + unique nFL reads\n")
	#f.write("# norm_nfl_amb: count_nfl_amb / total number of all reads\n")
	f.write("# Total Number of FL reads: {0}\n".format(total_reads))
	#f.write("# Total Number of FL + unique nFL reads: {0}\n".format(total_reads))
	#f.write("# Total Number of all reads: {0}\n".format(total_reads))
	f.write("#\n")
	f.write("pbid\tcount_fl\tnorm_fl\n") #f.write("pbid\tcount_fl\tcount_nfl\tcount_nfl_amb\tnorm_fl\tnorm_nfl\tnorm_nfl_amb\n")
	
def make_collapsed_abundance_txt(collapsed_group_txt, HQ_reads):
	dirc = os.path.dirname(collapsed_group_txt)
	basename = os.path.basename(collapsed_group_txt)
	newfile = open(dirc+'/'+basename.replace('group','abundance'),'w')
	
	try:
		HQ_reads = check_fasta.check_fasta_fmt(HQ_reads)
	except:
		if not os.path.exists(HQ_reads):
			sys.exit(-1)
			
	total_reads = int(os.popen('wc -l %s' %(HQ_reads)).read().split()[0])/2
	
	write_header(newfile, total_reads)
	
	a=b=c=d=e=f=g=h=j=k=l=n=m=x = 0
	
	file2 = open(collapsed_group_txt, 'r')
	length = int(os.popen('wc -l %s' %(collapsed_group_txt)).read().split()[0])
	print('\nWe are looping over ' + str(length) + ' groups\n')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x<length:
		group = file2.readline().strip()
		split_group = group.split()
		group_name = split_group[0]
		reads = split_group[1].split(',')
		newfile.write(group_name + '\t' + str(len(reads)) + '\t' + '0' + '\t' + '0' + '\t' + "{:.2e}".format(len(reads)/total_reads) + '\t' + '0' + '\t' + '0' + '\n')
		x += 1
		bar.update(x)
	bar.finish()
	newfile.close()
	print('\nDone, bye \n')
	
def make_collapsed_abundance_txt2(all_samples_chained_count_txt, total_reads):
	dirc = os.path.dirname(all_samples_chained_count_txt)
	basename = os.path.basename(all_samples_chained_count_txt)
	newfile = open(dirc+'/flair.firstpass.q.counts_NEW','w') #all_samples.chained.abundance.txt','w')
	
	write_header(newfile, total_reads)
	'''
	for line in open(all_samples_chained_count_txt, 'r'):
		if not line.startswith('superPBID'):
			items = line.split()
			items2 = list(filter(lambda a: a != 'NA', items[1:]))
			abundance = round(eval('+'.join(items2))/len(items2))
			#newfile.write(items[0] + '\t' + str(abundance) + '\t' + '0' + '\t' + '0' + '\t' + "{:.2e}".format(abundance/total_reads) + '\t' + '0' + '\t' + '0' + '\n')
			newfile.write(items[0] + '\t' + str(abundance) + '\t' + "{:.2e}".format(abundance/total_reads) + '\n')
	'''
	for line in open(all_samples_chained_count_txt, 'r'):
		if not line.startswith('superPBID'):
			ID = line.split()[0]
			count = line.strip().split()[1]
			newfile.write(ID + '\t' + str(round(float(count))) + '\t' + "{:.2e}".format(round(float(count))/total_reads) + '\n')
	newfile.close()
			
def get_reverse(fasta, fmt):
	dirc = os.path.dirname(fasta)
	basename = os.path.basename(fasta)
	newfile = open(dirc+'/'+'reverse_'+basename,'w')
	
	try:
		fasta = check_fasta.check_fasta_fmt(fasta)
	except:
		if not os.path.exists(fasta):
			sys.exit(-1)
			
	total_reads = int(os.popen('wc -l %s' %(fasta)).read().split()[0])/2
	print('We are looping over ' + str(total_reads) + ' reads\n')
	
	length=0
	x = 0
	for line in open(fasta):
		length+=1
		
	ofasta = open(fasta, 'r')
	
	a=b=c=d=e=f=g=h=j=k=l=n=m = 0
	complement = {'a':'T','c':'G','g':'C','t':'A','n':'N','A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length:
		seq=''
		rd_name = ofasta.readline()
		rd_seq = ofasta.readline().strip()
		if rd_name.strip().split('.')[1:][0] == '+f' or  rd_name.strip().split('.')[1:][1] == 'fb':
			newfile.write(rd_name+rd_seq+'\n')
			a += 1
		elif rd_name.strip().split('.')[1:][0] == '-r' or  rd_name.strip().split('.')[1:][2] == 'rb':
			for item in rd_seq[::-1]:
				seq=seq+complement[item]
			newfile.write(rd_name+seq+'\n')
			b += 1
		else:
			newfile.write(rd_name+rd_seq+'\n')
			c += 1
		x+=2
		bar.update(x)
	bar.finish()
	newfile.close()
	print('\nForward reads are ' + str(a))
	print('\nReverse reads are ' + str(b))
	print('\nReads with unknown direction are ' + str(c)+'\n')

def manip_gbk(gbk):
	dirc = os.path.dirname(gbk)
	basename = os.path.basename(gbk)
	newfile = open(dirc+'/'+basename+'.fasta','w')
	
	total_lines = int(os.popen('wc -l %s' %(gbk)).read().split()[0])
	print('We are looping over ' + str(total_lines) + ' lines\n')
	
	ofasta = open(gbk, 'r')
	x = 0
	xp = ''
	length = 0
	defn = ''
	gene = ''
	db_xref = ''
	seq1 = []
	seq = ''
	ref= ''
	ORIGIN = 0
	END = 0
	bar = progressbar.ProgressBar(maxval=total_lines, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < total_lines:
		#line = ofasta.readline()
		#print([line])
		line = ofasta.readline().strip()
		line2 = line.split()
		#if len(line2) > 0:# and line2[0] == '//':
		#	print([line2])
		if len(line2) >= 3 and line2[0] == 'LOCUS':
			xp = line2[1]
			length = line2[2]
		elif len(line2) >= 4 and line2[0] == 'DBSOURCE' and line2[3].startswith('NC_'):
			ref = line2[3].strip('()')
		elif len(line2) >= 1 and line2[0].startswith('/product='):
			#print(line+'\n'+line2[0])
			defn = [line][0].lstrip('/product=').strip('"')
		elif len(line2) >= 1 and line2[0].startswith('/gene='):
			gene = line2[0].strip('/gene=').strip('"')
		elif len(line2) >= 1 and line2[0].startswith('/db_xref="GeneID'):
			db_xref = line2[0].strip('/db_xref="GeneID:').strip('"')
		elif len(line2) >= 1 and line2[0].startswith('(NW_'):
			ref = line2[0].strip('()')
		elif len(line2) == 1 and line2[0] == 'ORIGIN' and ORIGIN == 0:
			#print(line2)
			ORIGIN = 1
		elif ORIGIN == 1 and END == 0 and line2[0] != '//':
			#print(seq1)
			#print([line][0])
			seq1 += [[line][0]]
		elif len(line2) > 0 and line2[0] == '//':
			#print(1)
			ORIGIN = 0
			for i in seq1:
				i = i.strip('1234567890//').replace(' ','')
				seq = seq+i
			newfile.write('|'.join(['>'+db_xref,ref,gene,length,xp,defn]) + '\n' + seq + '\n')
			#print('|'.join(['>'+db_xref,ref,gene,length,xp,defn]) + '\n' + str(seq) + '\n')
			xp = ''
			length = 0
			defn = ''
			gene = ''
			db_xref = ''
			seq1 = []
			seq = ''
			ref= ''
		x += 1
		bar.update(x)
	bar.finish()
	newfile.close()
	
def man_uniprot_dat(dat_file):
	dirc = os.path.dirname(dat_file)
	basename = os.path.basename(dat_file)
	newfile = open(dirc+'/'+basename+'filtered','w')
	
	total_lines = int(os.popen('wc -l %s' %(dat_file)).read().split()[0])
	print('\nWe are looping over ' + str(total_lines) + ' lines\n')
	
	ofasta = open(dat_file, 'r')
	x = 0
	ID = '.'
	new_protein = 'no'
	AC = '.'
	OC = []

	bar = progressbar.ProgressBar(maxval=total_lines, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < total_lines:
		line = ofasta.readline().strip()
		line2 = line.split()

		if line2[0] == 'ID':
			ID = line2[1]
		elif line2[0] == 'AD':
			AC = line2[1].strip(';')
		elif line2[0] == 'OC':
			OC += [''.join(line2[1:])]
		
		if line2[0] == "//":
			new_protein = 'yes'
			newfile.write('\t'.join([ID,AC,''.join(OC)])+'\n')
			ID = '.'
			AC = '.'
			OC = []
		else:
			new_protein = 'no'
		x += 1
		bar.update(x)
	bar.finish()
	newfile.close()
			
def make_fasta_single_line(fasta):
	dirc = os.path.dirname(fasta)
	basename = os.path.basename(fasta)
	newfile = open(dirc+'/'+'single_'+basename,'w')
	
	#try:
	#	fasta = check_fasta.check_fasta_fmt(fasta)
	#except:
	#	if not os.path.exists(fasta):
	#		sys.exit(-1)
			
	total_lines = int(os.popen('wc -l %s' %(fasta)).read().split()[0])
	print('We are looping over ' + str(total_lines) + ' lines\n')
	x = 0
	ofasta = open(fasta, 'r')
	bar = progressbar.ProgressBar(maxval=total_lines, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < total_lines:
		'''
		seq_id = ofasta.readline().strip()		
		seq = ofasta.readline()		
		newfile.write(seq_id + '\t' + seq)
		x += 2
		'''		
		seq_id = ofasta.readline()
		s_id = seq_id.split('\t')[0]
		seq = seq_id.split('\t')[-1]
		newfile.write(s_id + '\n' + seq)		
		x += 1
		#'''
		bar.update(x)
	bar.finish()
	newfile.close()
			
def add_uniprot_id_to_clusters(clusters, uniprot):
	uniprot_dict = {}
	
	x = 0
	
	length = int(os.popen('wc -l %s' %(uniprot)).read().split()[0]) #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
	
	ofile2=open(uniprot,'r')
	print('\nCreating dictionary for %s lines\n' %(str(length)))
	bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length: #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
		#'''
		a=ofile2.readline().strip()
		if a.startswith('gene'):
			d=a.split("|")[0]
			e=a.split("|")[2]
			uniprot_dict.setdefault(d.strip(), [a,e])
		elif a.startswith('PB'):
			a1 = a.split()
			uniprot_dict.setdefault(a1[0].strip(), [a,a])
		x+=1
		bar.update(x)
		#'''
	bar.finish()
	print(len(uniprot_dict.keys()))
	print('\nDone creating dictionary, moving on to querrying...\n')
	#'''
	r = 0
	print(len(uniprot_dict.keys()))
	for x,y in uniprot_dict.items():
		if r < 1: #if x.startswith('PB.9207.3') and r < 1:
			print([x],y)
			r +=1
	#'''
	
	dirc = os.path.dirname(clusters)
	name = os.path.basename(clusters)
	nfile1 = open(dirc+'/'+name.replace('.txt','_plus_dmel_homologue.txt'), 'w') #open(name.replace('_2','.gbk.fasta'), 'w')
	nfile2 = open(dirc+'/'+'clusters_without_dmel_homologue', 'w') #open(name.replace('_2','.gbk.nofasta'), 'w')
	
	a=b=c=d=e=f=g = 0
	
	length = int(os.popen('wc -l %s' %(clusters)).read().split()[0])
	print('We are querying dictionary for ' + str(length) + ' items\n')
	cluster_lines = open(clusters, 'r')
	#next(cluster_lines)
	cluster_lines = cluster_lines.readlines()
	cluster_lines = cluster_lines[1:]
	new_clusters = {}
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in cluster_lines:
		cluster_n = line.strip().split()[0]
		line1 = line.strip().split()[1].split('_')[0].strip('.path1')
		line2 = line.strip().split()[1].split('_')[1]
		#if line1 == 'PB.9207.3':
		#	print('yes')
		if line1 in uniprot_dict.keys():
			uniprot_ID1 = uniprot_dict[line1][0].split()[1].split('|')[1]
			uniprot_ID2 = uniprot_dict[line1][0].split()[1].split('|')[2]
			nfile1.write('\t'.join([line.strip(),uniprot_ID1,uniprot_ID2])+'\n')
				#del uniprot_dict[line1]		
			a += 1
			if cluster_n not in new_clusters.keys():
				new_clusters.setdefault(cluster_n,[uniprot_ID2.replace('_DROME','')])
			elif cluster_n in new_clusters.keys():
				new_clusters[cluster_n].append(uniprot_ID2.replace('_DROME',''))
		else:
			nfile2.write('\t'.join([line.strip(),line1,line2])+'\n')
			b += 1
		c+=1		
		bar.update(c)
	for x,y in new_clusters.items():
		cluster_file = open(dirc+'clusters/cluster_'+x, 'w')
		for z in y:
			cluster_file.write(z+'\n')
		cluster_file.close()
			
	bar.finish()
	nfile1.close()
	nfile2.close()
	
	print('\n'+str(a)+'\t'+str(b)+'\n')
	print('Found ' + str(a) + ' matches which is %s percent of all the queries u gave...' %(str("{0:.1f}".format(round(100*a/length),1))))
	
def manip_samfile(samfile):
	dirc = os.path.dirname(samfile)
	name = os.path.basename(samfile)
	a=b=c=d=e=f=g=x = 0
	sam_dic = {}
	#nfile1 = open(dirc+'/'+name+'_edited', 'w')
	nfile1 = open(dirc+'/'+name.replace('sam','fusions'), 'w')
	length = int(os.popen('wc -l %s' %(samfile)).read().split()[0])
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(samfile):
		alignment = line.split('\t')
		'''
		#if alignment[2].startswith('NW_') or alignment[2].startswith('NC_'):
		if alignment[1].split(":")[1].startswith('NW_') or alignment[1].split(":")[1].startswith('NC_'):
			nfile1.write(line.replace(alignment[1].split(":")[1],alignment[1].split(":")[1].replace('.','_')))
		else:
			nfile1.write(line)
		'''
		'''
		extra_fields = alignment[11:]
		aln_type = '*'
		for i in extra_fields:
			if len(extra_fields)>1 and 'tp:A:' in i:
				aln_type = i.strip('tp:A:')
				if aln_type == 'P':
					nfile1.write(line)
					a += 1
				elif aln_type == 'S':
					b += 1
				elif aln_type == 'I':
					c += 1
				else:
					e += 1
		'''
		#get fusion genes
		if line.startswith('@') or alignment[2] == '*':
			continue
		else:
			read = alignment[0]
			ref = alignment[2].split('_')[-1]
			if read not in sam_dic.keys():
				sam_dic.setdefault(read,[ref])
			elif read in sam_dic.keys():
				sam_dic[read].append(ref)
		d+=1		
		bar.update(d)
	bar.finish()
	for read,ref in sam_dic.items():
		ref = list(set(ref))
		if len(ref) < 2:
			continue
		elif len(list(set(ref))) == 2:
			nfile1.write('\t'.join([read, str(len(ref)), '--'.join(ref)])+'\n')
		elif len(ref) > 2:
			x = 0
			itern = len(ref)
			fusions = []
			while x < itern:
				curr_gene = ref[x]
				for i in ref[x+1:]:
					fusions += [curr_gene+'--'+i]
				x += 1
			for fusion in fusions:
				nfile1.write('\t'.join([read, 'x', fusion+'\n']))
				
	nfile1.close()
	print('Found ' + str(a) + ' primary alignments which is %s percent of all the lines u gave...' %(str("{0:.1f}".format(round(100*a/length),1))))
	print('Found ' + str(b) + ' secondary alignments')
	print('Found ' + str(c) + ' Inversions alignments')
	print('Found ' + str(e) + ' unknown alignments')
	print(length)

def manip_samfile2(samfile, details):
	dirc = '/home/banthony/scratch/analysis/covid19/B004/' #os.path.dirname(samfile)
	name = os.path.basename(samfile)
	a=b=c=d=e=f=g=x = 0
	sam_dic = {}
	det_dic = {}
	headers = []

	for line in open(details, 'r'):
		barcode = line.strip().split('\t')[-1].replace('BC','')
		subfolder = line.strip().split('\t')[0].split('_')[-2]
		if 'v3' in line:
			folder = 'artic_v3'
		elif 'v4' in line:
			folder = 'artic_v4'
		elif 'Midnight' in line:
			folder = 'midnight'
		elif 'SNAP' in line:
			folder = 'snap'
		det_dic.setdefault(barcode,folder+'/'+subfolder)
	
	for line in open(samfile, 'r'):
		if line.startswith('@'):
			headers += [line]
		barcode = line.split('_')[0].replace('barcode','')
		if barcode not in sam_dic.keys():
			sam_dic.setdefault(barcode, [line])
		else:
			sam_dic[barcode].append(line)
	
	#bar = progressbar.ProgressBar(maxval=len(det_dic.keys), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	#bar.start()
	for barcodex, folderx in det_dic.items():
		print([barcodex,folderx])
		name = 'barcode' + barcodex+'.sam'
		nfile1 = open(dirc + folderx + '/'+name, 'w')
		nfile1.write(''.join(headers))
		nfile1.write(''.join(sam_dic[barcodex]))
		nfile1.close()
		a += 1
		#bar.update(a)
	#bar.finish()
	
def add_stuff(ptn_len, homolog):
	
	dirc = os.path.dirname(homolog)
	name = os.path.basename(homolog)
	
	homolog = open(homolog,'r')
	homolog = homolog.readlines()
	
	nfile1 = open(dirc+'/'+name + '.new', 'w')
	nfile1.write('\t'.join(['qseqid','sseqid','pident','aln_length','qstart','qend','sstart','send','qlen','slen','qcov','scov','evalue'])+'\n')
	x=k = 0
	length = int(os.popen('wc -l %s' %(ptn_len)).read().split()[0])
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(ptn_len):
		pname = line.split('|')[1]+'|'+line.split('|')[2]
		for hom in homolog:
			if pname in hom:
				try:
					hom_spt = hom.split()
					#print(hom_spt[0])
					qlen = hom_spt[0].split('|')[4]
					slen = line.strip().split('|')[-1]
					perc_qcov = round(100*(int(hom_spt[3])-int(hom_spt[4])-int(hom_spt[5]))/int(qlen))
					perc_scov = round(100*(int(hom_spt[3])-int(hom_spt[4])-int(hom_spt[5]))/int(slen))
					nfile1.write('\t'.join([hom_spt[0],line.strip(),hom_spt[2],hom_spt[3],hom_spt[6],hom_spt[7],hom_spt[8],hom_spt[9],\
											str(qlen),str(slen),str(perc_qcov),str(perc_scov),hom_spt[10]])+'\n')
				except:
					k += 1
		x+=1
		bar.update(x)
	bar.finish()
	nfile1.close()
	print(k)
	
def man_gft(gtf, intron_counts):
	#This code takes a gtf and will return the transcripts associated with each gene such that: gene1 rna1,rna2,rna3
	dirc = os.path.dirname(gtf)
	name = os.path.basename(gtf)
	x=0
	gene_dic = {}
	rna_dic = {}
	gene_id = ''
	transcript_id = ''
	length = int(os.popen('wc -l %s' %(gtf)).read().split()[0])
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(gtf, 'r'):
		line = line.split('\t')
		line2 = line[-1].split(';')
		for i in line2:
			#print(i)
			if i.strip().startswith('transcript_id'):
				transcript_id = i.strip('transcript_id ')
				transcript_id = transcript_id.strip('"')
			elif i.strip().startswith('gene_id'):
				gene_id = i.strip('gene_id ')
				gene_id = gene_id.strip('"')

		if gene_id != '' and transcript_id != '' and gene_id not in gene_dic.keys():
			gene_dic.setdefault(gene_id,[transcript_id])
		elif gene_id != '' and transcript_id != '' and gene_id in gene_dic.keys():
			gene_dic[gene_id].append(transcript_id)
		elif gene_id == '' and transcript_id != '' and gene_id not in gene_dic.keys():
			gene_dic.setdefault(transcript_id,[transcript_id])
		elif gene_id == '' and transcript_id != '' and gene_id in gene_dic.keys():
			gene_dic[transcript_id].append(transcript_id)
		
		x+=1
		bar.update(x)
	bar.finish()
				
	nfile1 = open(dirc+'/'+name + '.new', 'w')
	for x,y in gene_dic.items():
		nfile1.write(x+'\t'+','.join(list(set(y)))+'\n')
		for i in list(set(y)):
			rna_dic.setdefault(i,x)
	nfile1.close()
	x=0
	#This part of the code is to add the gene id to the file that has transcript IDs
	dirc = os.path.dirname(intron_counts) 
	name = os.path.basename(intron_counts)
	nfile1 = open(dirc+'/'+name + '.new', 'w')
	length = int(os.popen('wc -l %s' %(intron_counts)).read().split()[0])
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(intron_counts, 'r'):		
		line2 = line.split('\t')
		if line2[3] in rna_dic.keys():
			nfile1.write(line.strip() + '\t' + rna_dic[line2[3]] + '\n')
		x+=1
		bar.update(x)
	bar.finish()
	nfile1.close()
	
def man_gft2(gtf):
	#This code takes a gtf and will return the transcripts that have mRNA and their associated info such that: transcript_id gene_id gene_symbol gene_length mRNA_length strand
	dirc = os.path.dirname(gtf)
	name = os.path.basename(gtf)
	nfile1 = open(dirc+'/'+name + '.eddited', 'w')
	a=b=c=d=e=f=0
	gene_dic = {}
	transcript_dic = {}
	gene_id = ''
	transcript_id = ''
	length = int(os.popen('wc -l %s' %(gtf)).read().split()[0])
	verified = ['PB.107','PB.157','PB.165','PB.166','PB.170','PB.174','PB.269','PB.360','PB.394','PB.417','PB.428','PB.439','PB.599','PB.616','PB.635','PB.65','PB.676']
	gene_IDPat = re.compile(r'(gene_id "[^"]*)"')
	trans_IDPat = re.compile(r'transcript_id "([^"]+)"')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(gtf, 'r'):
		#print(line.split('\t')[-1])
		line1 = line.split('\t')
		line2 = line[-1].split(';')
		sq_gene_id = gene_IDPat.search(line).group(1)
		trans_ID = trans_IDPat.search(line).group(1)
		tofu_gene_id = trans_ID.rstrip('.'+trans_ID.split('.')[-1])
		if tofu_gene_id in verified:
			to_write = line.replace(sq_gene_id,'gene_id "'+tofu_gene_id)
			nfile1.write(to_write.replace("PB.","PB.CAJ."))
		else:
			to_write = line.replace(sq_gene_id,'gene_id "'+trans_ID)
			nfile1.write(to_write.replace("PB.","PB.CAJ."))
		'''
		if line1[2] == "gene":
			#genlen = int(line1[4]) - int(line1[3])			
			#gene_IDPat = re.compile(r'gene_id "([^"]+)"')  #re.compile(r'gene_id "(.*?)"')
			#gene_symbPat = re.compile(r'gene_symbol "([^"]+)"')
			#sq_gene_id = gene_IDPat.search(line).group(1)
			#gene_symb = gene_symbPat.search(line).group(1)
			#trans_ID = trans_IDPat.search(line).group(1)
			#tofu_gene_id = trans_ID.rstrip('.'+trans_ID.split('.')[-1])
			if gene_id not in gene_dic.keys():
				gene_dic.setdefault(gene_id,[gene_id,gene_symb,str(genlen)])
			else:
				a += 1
		if line1[2] == "mRNA":
			mrnalen = int(line1[4]) - int(line1[3])
			mrnaDir = line1[6]
			gene_IDPat = re.compile(r'gene_id "([^"]+)"')
			transcript_IDPat = re.compile(r'transcript_id "([^"]+)"')
			gene_id = gene_IDPat.search(line).group(1)
			transcript_id = transcript_IDPat.search(line).group(1)
			if transcript_id not in transcript_dic.keys():
				transcript_dic.setdefault(transcript_id,gene_dic[gene_id]+[str(mrnalen),mrnaDir])
			else:
				b += 1
		'''
		c+=1
		bar.update(c)
	bar.finish()
	'''
	nfile1 = open(dirc+'/'+name + '.gene_mrna_info', 'w')
	for transcript,info in transcript_dic.items():
		nfile1.write(transcript+'\t'+'\t'.join(info)+'\n')
	'''
	nfile1.close()
	print(str(a)+'\t'+str(b))
	
def man_gft3(gtf, file1):
	#This code takes a gtf and will return the transcripts that have mRNA and their associated info such that: transcript_id gene_id gene_symbol gene_length mRNA_length strand
	dirc = os.path.dirname(gtf)
	name = os.path.basename(gtf)
	a=b=c=d=e=f=0
	gene_dic = {}
	transcript_dic = {}
	gene_id = transcript_id = ''	
	error_lines = []
	length = int(os.popen('wc -l %s' %(gtf)).read().split()[0])
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for linex in open(file1, 'r'):
		slinex = linex.split()
		pb_gene_id_items = slinex[0].split('.'); pb_gene_id = '.'.join(pb_gene_id_items[:2])
		if slinex[1].startswith('CAJ') and pb_gene_id not in gene_dic.keys():
			gene_dic.setdefault(pb_gene_id,[ [slinex[0]],[slinex[25]] ])
		elif slinex[1].startswith('CAJ') and pb_gene_id in gene_dic.keys():
			gene_dic[pb_gene_id][0].append(slinex[0])
			gene_dic[pb_gene_id][1].append(slinex[25])
		elif slinex[1].startswith('NW_') and slinex[7] == "novel" and pb_gene_id not in gene_dic.keys():
			gene_dic.setdefault(pb_gene_id,[ [slinex[0]],[slinex[25]] ])
		elif slinex[1].startswith('NW_') and slinex[7] == "novel" and pb_gene_id in gene_dic.keys():
			gene_dic[pb_gene_id][0].append(slinex[0])
			gene_dic[pb_gene_id][1].append(slinex[25])
		else:
			print(linex)
	#for gene_idx,itemx in gene_dic.items():
		
	for line in open(gtf, 'r'):
		line1 = line.split()
		line2 = line[-1].split(';')
		'''
		if line1[2] == "gene":
			genlen = int(line1[4]) - int(line1[3])			
			gene_IDPat = re.compile(r'gene_id "([^"]+)"')  #re.compile(r'gene_id "(.*?)"')
			gene_symbPat = re.compile(r'gene_name "([^"]+)"')
			gene_id = gene_IDPat.search(line).group(1).replace('path1','')
			gene_symb = gene_symbPat.search(line).group(1)
			if gene_id not in gene_dic.keys():
				gene_dic.setdefault(gene_id,[gene_id,gene_symb,str(genlen)])
			else:
				a += 1
		'''
		if line1[2] == "transcript":
			mrnalen = int(line1[4]) - int(line1[3])
			mrnaDir = line1[6]
			chrom = line1[0]
			gene_IDPat = re.compile(r'gene_id "([^"]+)"')
			transcript_IDPat = re.compile(r'transcript_id "([^"]+)"')
			gene_id = gene_IDPat.search(line).group(1).replace('mrna1','')
			transcript_id = transcript_IDPat.search(line).group(1)
			chrom_gene = chrom+'\t'+gene_id
			if chrom_gene not in gene_dic.keys():
				gene_dic.setdefault(chrom_gene,[[line1[3],line1[4]],[line1[6]]])
			else:
				#gene_dic[transcript_id][0].append( line1[3],line1[4] )
				gene_dic[chrom_gene][0]+= [ line1[3],line1[4] ]
				b += 1		
		c+=1
		bar.update(c)
	bar.finish()
				
	nfile1 = open(dirc+'/'+name.replace('gff','bed'), 'w')
	#nfile2 = open(dirc+'/'+name + '.errorLines', 'w')
	for transcript,info in gene_dic.items():
		nfile1.write(transcript.split()[0]+'\t'+ str(int(sorted(info[0])[0])-1000) + '\t' + str(int(sorted(info[0])[-1])+1000) + '\t' + transcript.split()[1]
					 + '\t' + '960' + '\t' + info[1][0] + '\n')
	#nfile2.write(''.join(error_lines))
	nfile1.close()
	#nfile2.close()
	print(str(a)+'\t'+str(b))
	#print(error_lines[1])
	
def get_genomic_ranges(file, readlength):
	dirc = os.path.dirname(file)
	name = os.path.basename(file)
	x=0
	gene_dic = {}
	length_dic = {}
	for i in open(readlength, 'r'):
		i = i.split()
		ref = i[6]
		length = i[8].strip("()bp. ")
		length_dic.setdefault(ref,length)
	length = int(os.popen('wc -l %s' %(file)).read().split()[0])
	#bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	#bar.start()
	for line in open(file, 'r'):
		line2 = line.split()
		key = ';'.join(line2[:2])
		value1 = line2[8]
		value2 = line2[9]
		if key not in gene_dic.keys():
			gene_dic.setdefault(key,[value1,value2])
		else:
			gene_dic[key].append(value1)
			gene_dic[key].append(value2)
	#print(gene_dic)
	nfile1 = open(dirc+'/'+name + '.new', 'w')
	for x,y in gene_dic.items():
		y = [int(x) for x in y]
		start = min(y)
		end = max(y)
		if int(start) >= 1000:
			start2 = str(int(start)-1000)
		else:
			start2 = '1'
		if x.split(';')[1] in length_dic.keys():
			end2 = length_dic[x.split(';')[1]]
		else:
			end2 = str(int(end)+1000)
		nfile1.write('\t'.join([ x.split(';')[1],start2,end2, ':'.join([x.split(';')[0],x.split(';')[1],str(start)+'_'+str(end)]), '.',start2,end2,x.split(';')[0], str((int(end)-int(start))) ])+'\n')
	nfile1.close()
	
def get_genomic_ranges2(file, gtf, ref):
	dirc = '/home/banthony/scratch/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/gfold/zygotic_genes_analysis/maternallydegradedCompletelyPromoter/' #os.path.dirname(file)
	name = os.path.basename(file)
	nfile1 = open(dirc+'/' + 'maternallydegradedCompletely.promoter.fasta', 'w')
	nfile2 = open(dirc+'/' + 'maternallydegradedCompletely.geneAndpromoter.fasta', 'w')
	nfile3 = open(dirc+'/' + 'maternallydegradedCompletely.geneSeq.fasta', 'w')
	a=b=c=d=e=f=g=h=i=j=k=l=n=m = 0
	gene_dic = {}
	gene_dic2 = {}
	gene_dic3 = {}
	ref_dic = {}
	for genex in open(gtf, 'r'):		
		#gene_IDPat = re.compile(r'gene_id "([^"]+)"')  #re.compile(r'gene_id "(.*?)"')
		#gene = gene_IDPat.search(genex).group(1)
		if len(genex.split()) > 5 and genex.split()[2] == 'gene':
			gene = genex.split()[-1].split(';')[0].strip('ID=')
		else:
			continue
		if gene not in gene_dic.keys():
			direction = genex.split()[6]
			gene_dic.setdefault(gene,direction)
			gene_dic2.setdefault(gene, [genex.split()[3],genex.split()[4]])
			gene_dic3.setdefault(gene,genex.split()[0])
			a += 1
	seqname = ''
	fasta = check_fasta.check_fasta_fmt(ref)
	for liney in open(fasta, 'r'):
		if liney.startswith('>'):
			seqname = liney.strip('>\n')
		elif not liney.startswith('>'):
			ref_dic.setdefault(seqname,liney.strip('\n'))
			seqname = 2
			b += 1
	for line in open(file, 'r'):
		line2 = line.split() #.split('_')
		gene = line2[0] #.replace('.path1','')
		try:
			direction = gene_dic[gene]
		except:
			print(gene + ' not in gtf')
			continue
			
		ref = line2[1].split(':')[0]
		ref = gene_dic3[gene]
		#coords = line2[1].split(':')[1].split('-')
		coords = gene_dic2[gene]
		min1 = int(coords[0])
		max1 = int(coords[1])
		genelength = (max1 - min1)+1		
		try:
			ref_dic[ref]
		except:
			print('Ref ' + ref + ' not found')
			continue
		if direction == '-':
			seq = ref_dic[ref]
			lenSeq = len(seq)
			
		if direction == '+' and min1 >= 1000:
			genePromoter = ref_dic[ref][min1-1000:min1]
			geneAndpromoter = ref_dic[ref][min1-1000:max1]
			geneSeq = ref_dic[ref][min1:max1]
			c += 1
		elif direction == '+' and min1 < 1000 and min1 > 7:
			genePromoter = ref_dic[ref][:min1]
			geneAndpromoter = ref_dic[ref][:max1]
			geneSeq = ref_dic[ref][min1:max1]
			d += 1
		elif direction == '+' and min1 < 1000 and min1 <= 7:
			genePromoter = ''
			geneAndpromoter = ''
			geneSeq = ref_dic[ref][min1:max1]
			e += 1
		elif direction == '-' and (lenSeq-max1) >= 1000:
			genePromoter = get_reverse_complement(ref_dic[ref][max1:(max1+1000)])
			geneAndpromoter = get_reverse_complement(ref_dic[ref][min1:(max1+1000)])
			geneSeq = get_reverse_complement(ref_dic[ref][min1:max1])
			f += 1
		elif direction == '-' and (lenSeq-max1) < 1000 and (lenSeq-max1) > 7:
			genePromoter = get_reverse_complement(ref_dic[ref][max1:])
			geneAndpromoter = get_reverse_complement(ref_dic[ref][min1:])
			geneSeq = get_reverse_complement(ref_dic[ref][min1:max1])
			g += 1
		elif direction == '-' and (lenSeq-max1) < 1000 and (lenSeq-max1) <= 7:
			genePromoter = ''
			geneAndpromoter = ''
			geneSeq = ref_dic[ref][min1:max1]
			h += 1
		else:
			print(gene + ' is weird\n')
			genePromoter = ''
			geneAndpromoter = ''
			geneSeq = ''
			i += 1
		genePromoter = genePromoter.strip('N')
		geneAndpromoter = geneAndpromoter.strip('N')
		geneSeq = geneSeq.strip('N')
		if len(genePromoter) > 7:
			nfile1.write(">"+gene+'\t'+ref+'\t'+direction+'\t'+'-'.join(coords)+'\n'+genePromoter+'\n') #nfile1.write(">"+gene+'\t'+ref+'\t'+direction+'\t'+line2[1].split(':')[1]+'\n'+genePromoter+'\n')
		if len(geneAndpromoter) > 25:
			nfile2.write(">"+gene+'\t'+ref+'\t'+direction+'\t'+'-'.join(coords)+'\n'+geneAndpromoter+'\n') #nfile2.write(">"+gene+'\t'+ref+'\t'+direction+'\t'+line2[1].split(':')[1]+'\n'+geneAndpromoter+'\n')
		if len(geneSeq) > 25:
			nfile3.write(">"+gene+'\t'+ref+'\t'+direction+'\t'+'-'.join(coords)+'\n'+geneSeq+'\n') #nfile3.write(">"+gene+'\t'+ref+'\t'+direction+'\t'+line2[1].split(':')[1]+'\n'+geneSeq+'\n')
	nfile1.close()
	nfile2.close()
	nfile3.close()
	print('\t'.join([str(a),str(b),str(c),str(d),str(e),str(f),str(g),str(h),str(i)]))
	
def pick_random(file):
	from random import randint
	dirc = '/home/banthony/scratch/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/gfold/zygotic_genes_analysis/' #os.path.dirname(file)
	name = os.path.basename(file)
	nfile1 = open(dirc+'/' + 'random.maternal_degraded_zygotic', 'w')
	gene_dic = {}
	for line in open(file, 'r'):
		gene = line.split()[0]
		if gene not in gene_dic.keys():
			gene_dic.setdefault(gene,line)
			
	random_n = randint(0,5)
	a=b = 0
	for gene,line in gene_dic.items():
		if a == random_n and b < 200:
			nfile1.write(line)
			a = 0
			b += 1
		else:
			a += 1
	nfile1.close()
	
def combine_lines(file1):
	dirc = os.path.dirname(file1)
	name = os.path.basename(file1)
	x=0
	reads_dic = {}
	for line in open(file1, 'r'):
		line2 = line.strip().split()
		key = line2[0]
		value1 = line2[1]
		if key not in reads_dic.keys():
			reads_dic.setdefault(key,[value1])
		else:
			reads_dic[key].append(value1)
	nfile1 = open(dirc+'/'+name + '.new', 'w')
	for x,y in reads_dic.items():
		nfile1.write('\t'.join([x,','.join(y)])+'\n')
	nfile1.close()

def bed12togtf(bed12):
	dirc = os.path.dirname(bed12)
	name = os.path.basename(bed12)
	nfile1 = open(dirc+'/'+name+'.new.gtf', 'w')
	#nfile2 = open(dirc+'/'+'gtf_new2', 'w')
	a=b=c=d=e=f=g= 0
	h=1
	length = int(os.popen('wc -l %s' %(bed12)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for i in open(bed12, 'r'):
		ftrs = i.strip().split('\t')
		left=ftrs[1]
		size=ftrs[10].split(',')
		start=ftrs[-1].split(',')
		for x in range(len(size)-1):
			#print(int(left))
			nfile1.write('\t'.join([ftrs[0],'Gnomon','exon',str(int(left)+int(start[x])+1),str(int(left)+int(start[x])+int(size[x])),'.',ftrs[5],'.','transcript_id "rna%s"; gene_id "%s";' %(str(h),ftrs[3])])+'\n')
			h+=1
		d+=1
		bar.update(d)
	bar.finish()
	nfile1.close()
	
def find_man_identifier(file):
	dirc = os.path.dirname(file)
	name = os.path.basename(file)
	nfile1 = open(dirc+'/'+name.replace('fa','ed.fa'), 'w')
	a=b=c=d=e=f=g=h= 0
	length = int(os.popen('wc -l %s' %(file)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(file, 'r'):
		locRegex = re.compile(r'\(LOC(\d+)\),')
		matches = locRegex.search(line)
		if matches:
			nfile1.write('>LOC'+matches.group(1)+'\n')
			h+=1
		else:
			nfile1.write(line)
		d+=1
		bar.update(d)
	bar.finish()
	nfile1.close()
	print('\nDiscovered ' + str(h) + ' lines' + ' which is %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*h/length),1))))
	
def alignqc_bed(file1):
	dirc = os.path.dirname(file1)
	name = os.path.basename(file1)
	x=y=0
	reads_dic = {}
	readsFreq_dic = {}
	merged = 'No'
	for line in open(file1, 'r'):
		line2 = line.strip().split()
		key = line2[3]
		value1 = line2[0].split("_")[-1] #value1 = line2[0]
		if key not in reads_dic.keys():
			reads_dic.setdefault(key,[value1])
		elif key in reads_dic.keys():
			reads_dic[key].append(value1)
	nfile1 = open(dirc+'/'+name + '.new', 'w')
	'''
	for read,genes in reads_dic.items():
		nfile1.write('\t'.join([read,'--'.join(genes)])+'\n')
	'''
	for read,genes in reads_dic.items():
		key2 = '--'.join(genes)
		if key2 not in readsFreq_dic.keys():
			readsFreq_dic.setdefault(key2,[read])
		elif key2 in readsFreq_dic.keys():
			readsFreq_dic[key2].append(read)
	for read,genes in readsFreq_dic.items():
		try:
			if read.split('--')[0] != read.split('--')[1]:
				nfile1.write('\t'.join([read,str(len(genes)),','.join(genes)])+'\n')
		except:
			print(read)
	
		'''
		y='--'.join(y)
		if "~" in y:
			merged = 'yes'
			one = y.split('--')
			if "~" in one[0] and "~" not in one[1] :
				sub = one[1]
				que = one[0]
			elif "~" not in one[0] and "~" in one[1] :
				sub = one[0]
				que = one[1]
			elif "~" in one[0] and "~" in one[1]:
				nfile1.write('\t'.join([x,y,merged])+'\n')
			print(y)
			print(que)
			que = que.split("~")
			for i in que:
				nfile1.write('\t'.join([x,'--'.join([i,sub]),merged])+'\n')
		elif "~" not in y:
			nfile1.write('\t'.join([x,y,merged])+'\n')
		'''
	nfile1.close()
	
def file_manip_u_fool(file1):
	dirc = os.path.dirname(file1)
	name = os.path.basename(file1)
	nfile1 = open(dirc+'/'+name + '.new', 'w')
	for line in open(file1, 'r'):
		if line.startswith(">"):
			nfile1.write("_".join(line.strip().split())+'\n')
		elif not line.startswith(">"):
			line = line.strip().strip("Nn")+'\n'
			nfile1.write(line)
	nfile1.close()

def collapse_overlaps(file1):
	newfile = open(file1+'.final_boleae_miRNAs.bed','w')	
	newfile.write('track name=boleae_miRNA description="miRNAs identified in B. oleae" useScore=1\n')
	file1_dic = {}
	w=x=y=z = 0
	ofile1 = open(file1,'r')
	rfile1 = ofile1.readlines()
	coord1x = 0
	coord2x = 0
	miRNAs = []
	miRNAs2 = []
	for line in range(len(rfile1)):
		#print([line,len(rfile1)])
		sline = rfile1[line].strip().split('\t')
		miRNA = '-'.join(sline[0].split("_")[0].split('-')[1:])
		coord1 = min(int(sline[5]), int(sline[6]))
		coord2 = max(int(sline[5]), int(sline[6]))
		contig = sline[1]
		if line <= (len(rfile1)-2):
			#print([line,len(rfile1)])
			#print(rfile1[line])
			sline_next = rfile1[line+1].strip().split('\t')
			miRNAn = '-'.join(sline_next[0].split("_")[0].split('-')[1:])
			coord1n = min(int(sline_next[5]), int(sline_next[6]))
			coord2n = max(int(sline_next[5]), int(sline_next[6]))
			contig2 = sline_next[1]
		#if miRNAn == '':
		#	print(sline_next[0])
		if coord1n-coord1 <=20 and coord2n-coord2 <=20 and line != (len(rfile1)-2) and contig == contig2:
			coord1x = min(coord1,coord1n)
			coord2x = max(coord2,coord2n)
			miRNAs += [miRNA,miRNAn]
			miRNAs2 += [sline[0],sline_next[0]]
		elif coord1n-coord1 <=20 and coord2n-coord2 <=20 and line == (len(rfile1)-2) and contig == contig2:
			coord1x = min(coord1,coord1n)
			coord2x = max(coord2,coord2n)
			miRNAs += [miRNA,miRNAn]
			miRNAs2 += [sline[0],sline_next[0]]
			miRNA_last = '_'.join(list(dict.fromkeys(miRNAs)))
			miRNAs2_last = ','.join(list(dict.fromkeys(miRNAs2)))
			if miRNA_last not in file1_dic.keys():
				file1_dic.setdefault(miRNA_last,['\t'.join([contig,str(coord1x),str(coord2x)]),miRNAs2_last])
			else:
				print(miRNA_last)
		elif coord1n-coord1 > 20 or coord2n-coord2 > 20 or line == (len(rfile1)-2) or contig != contig2:
			#miRNAs += [miRNA,miRNAn]
			#miRNAs2 += [sline[0],sline_next[0]]
			miRNA_last = '_'.join(list(dict.fromkeys(miRNAs)))
			miRNAs2_last = ','.join(list(dict.fromkeys(miRNAs2)))
			if miRNA_last == '':
				miRNA_last = miRNA
				miRNAs2_last = sline[0]
				coord1x = coord1
				coord2x = coord2
				#print(rfile1[line])
			if miRNA_last not in file1_dic.keys():
				file1_dic.setdefault(miRNA_last,['\t'.join([contig,str(coord1x),str(coord2x)]),miRNAs2_last])
			else:
				print(miRNA_last) #n = 0 #print(miRNA_last)
			coord1x = 0
			coord2x = 0
			miRNAs = []
			miRNAs2 = []
		else:
			print(rfile1[line])
	for mirna,details in file1_dic.items(): #make bedfile using format shown at https://genome.ucsc.edu/FAQ/FAQformat.html#format1
		#print([details[0]])
		coord1 = details[0].split('\t')[1]
		coord2 = details[0].split('\t')[2]
		newfile.write('\t'.join([details[0],mirna,'960','+',coord1,coord2,'255,0,0','1',str(int(coord2)-int(coord1)),'0',details[1]])+'\n')
		#newfile.write('\t'.join([mirna,details[0],details[1]])+'\n')
	newfile.close()
	#print('done')
	
def cant_use_this_again(file1):
	import numpy as np
	dirc = os.path.dirname(file1)
	name = os.path.basename(file1)
	newfile = open(dirc+'/'+name + '.dmelOrthos','w')
	a=0
	dic = {'bdo':[],'bol':[],'cca':[],'dme':[],'mdo':[],'zcu':[]}
	#for line in next(open(file1, 'r')):
	line = open(file1, 'r')
	next(line)
	for linen in line:
		a+=1
		#print(line)
		line = linen.split()
		
		bdo = int(line[1])
		bol = int(line[2])
		cca = int(line[3])
		dme = int(line[4])
		mdo = int(line[5])
		zcu = int(line[6])
		if bdo != 0 and bol != 0 and cca != 0 and dme != 0 and mdo != 0 and zcu != 0:
			dic['bdo']+=list(np.repeat(np.array(["A"+str(a)]), [bdo]))
			dic['bol']+=list(np.repeat(np.array(["A"+str(a)]), [bol]))
			dic['cca']+=list(np.repeat(np.array(["A"+str(a)]), [cca]))
			dic['dme']+=list(np.repeat(np.array(["A"+str(a)]), [dme]))
			dic['mdo']+=list(np.repeat(np.array(["A"+str(a)]), [mdo]))
			dic['zcu']+=list(np.repeat(np.array(["A"+str(a)]), [zcu]))
		else:
			dic['bdo']+=list(np.repeat(np.array(["A"+str(a)]), [bdo]))
			dic['bol']+=list(np.repeat(np.array(["A"+str(a)]), [bol]))
			dic['cca']+=list(np.repeat(np.array(["A"+str(a)]), [cca]))
			dic['dme']+=list(np.repeat(np.array(["A"+str(a)]), [dme]))
			dic['mdo']+=list(np.repeat(np.array(["A"+str(a)]), [mdo]))
			dic['zcu']+=list(np.repeat(np.array(["A"+str(a)]), [zcu]))
	bdo = open(file1+'bdo','w')
	bol = open(file1+'bol','w')
	cca = open(file1+'cca','w')
	dme = open(file1+'dme','w')
	mdo = open(file1+'mdo','w')
	zcu = open(file1+'zcu','w')
	for sps,ortho in dic.items():
		if sps == 'bdo':
			bdo.write('\n'.join(ortho))
		elif sps == 'bol':
			bol.write('\n'.join(ortho))
		elif sps == 'cca':
			cca.write('\n'.join(ortho))
		elif sps == 'dme':
			dme.write('\n'.join(ortho))
		elif sps == 'mdo':
			mdo.write('\n'.join(ortho))
		elif sps == 'zcu':
			zcu.write('\n'.join(ortho))
def inutile_xyz(gff3):
	import numpy as np
	dirc = os.path.dirname(gff3)
	name = os.path.basename(gff3)
	newfile1 = open(name + '.CountSummary','w')
	newfile2 = open(name + '.LengthSummary','w')
	for_counts = {}
	for_len = {}
	for line in open(gff3,'r'):
		if not line.startswith('##'): # and not 'match_part' in line:
			ID = line.split('\t')[-1].split('_')[0].strip('ID=') #.split('-')[0]
			match = line.split()[2]
			ref = line.split()[0]
			start = line.split()[3]
			end = line.split()[4]
			length = str((max(int(start),int(end)) - min(int(start),int(end)))+1)+'\n'
			strand = line.split()[6]
			details = '\t'.join(line.split('\t')[-1].split()[0].replace('Target=','').split(';')[-1].split(':'))
			TE = details.split('\t')[0]
			details = details+'\t'+'\t'.join(list(np.repeat("UN", (5-len(details.split())) )))
			details = '\t'.join([ref,start,end,strand,ID,details.strip('\t'),length])
			
			if '-' not in ID:
				for_counts.setdefault(ID,details)
				for_len.setdefault(ID,details)
			elif match == 'match_part' and '-' in ID:
				if ID.split('-')[0].replace('p','s') in for_len.keys():
					del for_len[ID.split('-')[0].replace('p','s')]
					for_len.setdefault(ID,details)
				else:
					for_len.setdefault(ID,details)
	for key,value in for_counts.items():
		if 'MITE-Hunter' in value:
			value1 = value.strip().split()
			newfile1.write('\t'.join(value1[:6])+'\t'+'\t'.join(['Class_II','MITE','NA','UN'])+'\t'+value1[-1]+'\n')
		elif 'SINE-finder' in value:
			value1 = value.strip().split()
			newfile1.write('\t'.join(value1[:6])+'\t'+'\t'.join(['Class_I','SINE','incomplete','UN'])+'\t'+value1[-1]+'\n')
		else:
			newfile1.write(value)
	for key,value in for_len.items():
		if 'MITE-Hunter' in value:
			value1 = value.strip().split()
			newfile2.write('\t'.join(value1[:6])+'\t'+'\t'.join(['Class_II','MITE','NA','UN'])+'\t'+value1[-1]+'\n')
		elif 'SINE-finder' in value:
			value1 = value.strip().split()
			newfile2.write('\t'.join(value1[:6])+'\t'+'\t'.join(['Class_I','SINE','incomplete','UN'])+'\t'+value1[-1]+'\n')
		else:
			newfile2.write(value)
	newfile1.close()
	newfile2.close()
	counts = {}
	length = {}
def nevertouseagain(fasta):
	dirc = os.path.dirname(fasta)
	name = os.path.basename(fasta)
	newfile1 = open(name + '.CountSummary','w')
	'''
	for line in open(fasta, 'r'):
		if line.startswith(">"):
			line1 = line.strip().split("/")
			if '_' in line1[2]:
				line2 = line1[2].split('_')[1]
			else:
				line2 = line1[2]
			if line2 == 'unknown' or line2 == 'NA' or line2 == 'incomplete':
				line2 = 'unknown'
			newfile1.write('/'.join([line1[0],line1[1],line2])+'\n')
		else:
			newfile1.write(line.replace('/',''))
	newfile1.close()
	'''
	for line in open(fasta, 'r'):
		if line.startswith(">"):
			line1 = line.strip().split(":")
			#MITE and SINE specific tools
			if 'MITE-Hunter' in line1[0]:
				newfile1.write(line1[0]+'#'+'DNA'+'\n')
			elif 'SINE-finder' in line1[0]:
				newfile1.write(line1[0]+'#'+'SINE'+'\n')
				#un classed
			elif 'Class_noCat' in line1[1] or 'Class_NA' in line1[1]:
				newfile1.write(line1[0]+'#'+'unknown'+'\n')
				#DNA transposons
			elif line1[2].split('|')[0] == 'Helitron' or line1[2].split('|')[0] == 'MITE' or line1[2].split('|')[0] == 'Maverick' or line1[2].split('|')[0] == 'TIR':
				newfile1.write(line1[0]+'#'+'DNA'+'\n')
			elif line1[1] == 'Class_II' and line1[2] == 'noCat':
				newfile1.write(line1[0]+'#'+'DNA'+'\n')
				#retrotransposons
			elif line1[2].split('|')[0] == 'DIRS' or line1[2].split('|')[0] == 'LARD' or line1[2].split('|')[0] == 'PLE' or line1[2].split('|')[0] == 'TRIM':
				newfile1.write(line1[0]+'#'+'LTR'+'\n')
			elif line1[1] == 'Class_I' and line1[2] == 'noCat':
				newfile1.write(line1[0]+'#'+'RNA'+'\n')
			elif line1[2].split('|')[0] == 'LINE' or line1[2].split('|')[0] == 'LTR' or line1[2].split('|')[0] == 'SINE':
				newfile1.write(line1[0]+'#'+line1[2].split('|')[0]+'\n')
			else:
				print(line)
		else:
			newfile1.write(line.replace('/',''))
	newfile1.close()
	
def	find_low_trembl(file1):
	dirc = os.path.dirname(file1)
	name = os.path.basename(file1)
	newfile1 = open('tr_and_sp_blastp_hit_comparison','w')
	gene_dic = {}
	a=b=c=d=e=t=s=0
	for line in open(file1, 'r'):
		sline = line.split()
		genename = sline[0].split('|')[0]
		if genename not in gene_dic.keys():
			t=s=0
		if sline[1].startswith('tr') and genename not in gene_dic.keys() and t==0:
			gene_dic.setdefault(genename,['tr'+'_'+sline[10]])
			t+=1
		elif sline[1].startswith('sp') and genename not in gene_dic.keys() and s==0:
			gene_dic.setdefault(genename,['sp'+'_'+sline[10]])
			s+=1
		elif sline[1].startswith('sp') and genename in gene_dic.keys() and s==0:
			gene_dic[genename].append('sp'+'_'+sline[10])
			s+=1
		elif sline[1].startswith('tr') and genename in gene_dic.keys() and t==0:
			gene_dic[genename].append('tr'+'_'+sline[10])
			t+=1
	print(list(gene_dic.items())[0])
	for genename,blast in gene_dic.items():
		tr_ev = 100
		sp_ev = 100
		for x in blast:
			if x.startswith('tr'):
				tr_ev = x.split('_')[-1]		
			elif x.startswith('sp'):
				sp_ev = x.split('_')[-1]
		if float(tr_ev) < float(sp_ev):
			newfile1.write('\t'.join([genename,'\t'.join(blast),'yes\n']))
		elif float(tr_ev) > float(sp_ev):
			newfile1.write('\t'.join([genename,'\t'.join(blast),'no\n']))
	newfile1.close()
	
def reverse_some_seqs(fasta):
	print('\nWorking on ' + fasta)
	dirc = os.path.dirname(fasta)
	name = os.path.basename(fasta); end = ''
	if name.endswith('.fasta') or name.endswith('.fa'):
		end = '.'+name.split('.')[-1]
		name = name.rstrip('.fasta')
	newfile = open(dirc+'/'+name+'.reverseCompliment'+end,'w')
	for line in open(fasta, 'r'):
		if line.startswith(">"):
			newfile.write(line)
		else:
			newfile.write(get_reverse_complement(line.strip())+'\n')
	newfile.close()

def manipulate_sam_lines(sam, ref):
	dirc = os.path.dirname(sam)
	name = os.path.basename(sam)
	newfile = open(dirc+'/'+'class_reads_edited.fasta','w')
	ofasta = open(ref, 'r')
	rfasta = ofasta.readlines()
	#print(rfasta)
	for line in open(sam, 'r'):
		#print(line+'\n\n\n')
		if not line.startswith('@') and line.split()[2] != "*":
			items = line.strip().split('\t')
			refStart = (int(items[3])-1)
			cigar = items[5]
			seq = items[9]
			#soft clipping
			spat = re.compile(r'(\d+)S', re.I)
			soft = spat.findall(items[5])
			x = lclip = rclip = 0
			if len(soft) == 2:
				lclip = int(soft[0])
				rclip = int(soft[1])
			else:
				x+=1
			#print(seq)
			aligned_seq = seq[lclip:(len(seq)-rclip)]
			#print(aligned_seq) #print([str(lclip)+':'+str(len(seq)-rclip)]) #rfasta[1][:refStart])
			new_seq = rfasta[1][:refStart] + aligned_seq
			newfile.write('>'+items[0]+'\n'+new_seq+'\n')
		else:
			print(line.split()[0])
	newfile.close()
	print(x)

def get_longest_orfx(orfFasta):
	dirc = os.path.dirname(orfFasta)
	name = os.path.basename(orfFasta)
	newfile = open(dirc+'/'+'longest_orf.fasta','w')
	dic = {}
	x=0
	length = int(os.popen('wc -l %s' %(orfFasta)).read().split()[0])
	file1=open(orfFasta,'r')
	while x < length:
		name=file1.readline()
		sequence=file1.readline()
		name = '>'+name.split(':')[1]+'\n'
		if not name in dic.keys():
			dic.setdefault(name,[sequence.strip()+'\n'])
		else:
			dic[name].append(sequence.strip()+'\n')
		x+=2
	
	for name,sequence in dic.items():
		newfile.write(name+sorted(sequence, key=len)[-1])
	newfile.close()
	
def see_how_bad(file1):
	dicx = {}
	for line in open(file1, 'r'):
		if line.split()[1] not in dicx.keys():
			dicx.setdefault(line.split()[1],[line.split()[2]])
		else:
			dicx[line.split()[1]].append(line.split()[2])
	for x,y in dicx.items():
		#if x == 'gene7923':
		#	print(y)
		if '+' in y and '-' in y:
			print(x)
def delete_now(filex):
	for line in open(filex, 'r'):
		gene_IDPat = re.compile(r'gene_id "([^"]+)"')
		trans_IDPat = re.compile(r'transcript_id "([^"]+)"')
		gene_ID = gene_IDPat.search(line).group(1)	   
		trans_ID = trans_IDPat.search(line).group(1)
		print(gene_ID); print(trans_ID)
	
def sqanti_gff_edit(sqanti_gtf, ncbi_gtf, classn, grouping):
	dirc = os.path.dirname(sqanti_gtf)
	newfile = open(sqanti_gtf+'_updated', 'w')
	gene_IDPat = re.compile(r'gene_id "([^"]+)"')
	trans_IDPat = re.compile(r'transcript_id "([^"]+)"') 
	ncbi_gtf_dic = {}
	ncbi_trans_dic = {}
	ncbi_CDS_dic = {}
	ncbi_gene_trans_dic = {}
	a=b=c=d=e=f=g=h=i =0
	for line in open(ncbi_gtf, 'r'):
		if len(line.split()) == 10:
			gene_ID = trans_ID = trans_IDPat.search(line).group(1)
		else:
			gene_ID = gene_IDPat.search(line).group(1)	   
			trans_ID = trans_IDPat.search(line).group(1)
		exon_start = min( line.split()[3],line.split()[4] )
		exon_end = max( line.split()[3],line.split()[4] )
		if gene_ID not in ncbi_gene_trans_dic.keys() and line.split()[2] == 'exon':
			ncbi_gene_trans_dic.setdefault(gene_ID,[trans_ID])
		elif gene_ID in ncbi_gene_trans_dic.keys() and line.split()[2] == 'exon' and trans_ID not in ncbi_gene_trans_dic[gene_ID]:
			ncbi_gene_trans_dic[gene_ID].append(trans_ID)
		if gene_ID not in ncbi_gtf_dic.keys():
			ncbi_gtf_dic.setdefault(gene_ID,[line.split()[6],line.split()[2]])
		elif gene_ID in ncbi_gtf_dic.keys() and line.split()[2] not in ncbi_gtf_dic[gene_ID]:
			ncbi_gtf_dic[gene_ID].append(line.split()[2])
		if trans_ID not in ncbi_trans_dic.keys() and line.split()[2] == 'exon':
			ncbi_trans_dic.setdefault(trans_ID,[exon_start+'_'+exon_end])
		elif trans_ID in ncbi_trans_dic.keys() and line.split()[2] == 'exon':
			ncbi_trans_dic[trans_ID].append(exon_start+'_'+exon_end)
		elif trans_ID not in ncbi_CDS_dic.keys() and line.split()[2] == 'CDS':
			ncbi_CDS_dic.setdefault(trans_ID,[exon_start+'_'+exon_end])
		elif trans_ID in ncbi_CDS_dic.keys() and line.split()[2] == 'CDS':
			ncbi_CDS_dic[trans_ID].append(exon_start+'_'+exon_end)
	grouping_dic = {}
	for line in open(grouping, 'r'):
		flair_scaf = ''
		oldID_list = ['NW_'+x if not x.startswith('gene') else x for x in line.split()[1].strip().split(',') ]
		if len(oldID_list) == 1:
			oldID = oldID_list[0]
			flair_scaf = oldID.split('_')[-1]
		else:
			for i in oldID_list:
				oldID = oldID_list[0]
				flair_scaf = oldID.split('_')[-1]
				if flair_scaf.startswith('gene'):
					break			
		grouping_dic.setdefault(line.split()[0],flair_scaf)
	classn_dic = {}
	for line in open(classn, 'r'):
		items = line.split()
		if items[0] not in classn_dic.keys() and a > 0:
			newID = items[0]; oldID = items[0]; flair_scaf = grouping_dic[newID.replace('_dup2','')]; coding = items[29]; exons = items[4]; structure = items[5]; sqanti_assoc_gene = items[6]; sqanti_assoc_trans = items[7]; strand = items[2]
			classn_dic.setdefault(newID,[ oldID, flair_scaf, coding, exons, structure, sqanti_assoc_gene, sqanti_assoc_trans, strand ])
		a += 1
	
	sqanti_gene_dic = {}
	sqanti_trans_dic = {}
	sqanti_transCDS_dic = {}
	sqanti_lines = {}
	for line in open(sqanti_gtf,'r'):
		gene_IDPat = re.compile(r'gene_id "([^"]+)"')
		trans_ID = trans_IDPat.search(line).group(1)
		tofu_gene_ID = '.'.join(trans_ID.split('.')[:2]) #trans_ID.replace('.'+trans_ID.split('.')[-1],'')
		exon_start = min( line.split()[3],line.split()[4] )
		exon_end = max( line.split()[3],line.split()[4] )
		if trans_ID == 'PB.9001.17_dup2' or trans_ID == 'PB.7374.1_dup2':
			continue
		#add all gff lines to dic
		if trans_ID not in sqanti_lines.keys():
			sqanti_lines.setdefault(trans_ID,[line])
		elif trans_ID in sqanti_lines.keys():
			sqanti_lines[trans_ID].append(line)
		#add all tofu genes and their associated transcripts according to tofu
		if tofu_gene_ID not in sqanti_gene_dic.keys():
			sqanti_gene_dic.setdefault(tofu_gene_ID,[trans_ID])
		elif tofu_gene_ID in sqanti_gene_dic.keys() and trans_ID not in sqanti_gene_dic[tofu_gene_ID]:
			sqanti_gene_dic[tofu_gene_ID].append(trans_ID)
		#add all tofu transcripts and their associated exon coordinates to dic
		if trans_ID not in sqanti_trans_dic.keys() and line.split()[2] == 'exon':
			sqanti_trans_dic.setdefault(trans_ID,[exon_start+'_'+exon_end])
		elif trans_ID in sqanti_trans_dic.keys() and line.split()[2] == 'exon':
			sqanti_trans_dic[trans_ID].append(exon_start+'_'+exon_end)
		#add tofu transcripts and their features whether transcript,exon,CDS
		if trans_ID not in sqanti_transCDS_dic.keys():
			sqanti_transCDS_dic.setdefault(trans_ID,[line.split()[2]])
		elif trans_ID in sqanti_transCDS_dic.keys() and line.split()[2] not in sqanti_transCDS_dic[trans_ID]:
			sqanti_transCDS_dic[trans_ID].append(line.split()[2])
	sub_sqanti_trans_dic = sqanti_trans_dic	
	verified_dic = {'PB.6683.2':'gene10484','PB.783.1':'gene1193','PB.1318.2':'gene2224','B.2725.1':'gene4378','PB.3480.3':'gene5529',
					'PB.4624.1':'gene7445','PB.5033.1':'gene8008','PB.6683.1':'gene10487','PB.6683.3':'gene10486','PB.8347.2':'gene13178','PB.8649.1':'gene13629',
					'PB.8906.1':'gene14081','PB.2231.2':'gene3665','PB.2451.1':'gene4018','PB.2571.1':'PB.2571','PB.2571.2':'PB.2571','PB.2614.5':'gene4190',
					'PB.2692.1':'PB.2692.1','PB.2705.1':'PB.2705.1','PB.3182.1':'PB.3182.1','PB.4152.1':'PB.4152.1','PB.4220.1':'gene6755','PB.5027.1':'gene8004',
					'PB.5330.1':'gene8536','PB.5380.2':'PB.5380.2','PB.5383.1':'gene8614','PB.6019.1':'gene9442','PB.6140.2':'gene9693','PB.6435.1':'PB.6435.1',
					'PB.6464.1':'PB.6464.1','PB.6502.1':'gene10221','PB.8149.2':'gene12914','PB.8965.1':'gene14177'}
	genic_verified = {'PB.931.1':'gene1483','PB.1389.1':'gene2314','PB.1627.1':'gene2684','PB.2151.1':'gene3499','PB.3463.2':'gene5495','PB.4417.2':'gene7084',
					  'PB.4473.2':'gene7173','PB.4505.1':'gene7256','PB.4567.1':'gene7324','PB.4977.3':'gene7958','PB.5254.2':'gene8389','PB.5829.1':'gene9184',
					  'PB.5880.1':'gene9253','PB.6450.2':'gene10165','PB.8206.1':'gene2970'}
	verified_dic = genic_verified = {'PB.3742.1':'gene3018','PB.3742.2':'gene3018'}
	all_overs = {} #; print(sqanti_gene_dic); print('No no no')
	significant_overs = {}
	
	for newID,items in classn_dic.items():
		tofu_gene_IDx = '.'.join(newID.split('.')[:2]) #newID.replace('.'+newID.split('.')[-1],'')
		oldID = items[0]; flair_scaf = items[1]; coding = items[2]; exons = int(items[3]); structure = items[4]; sqanti_assoc_gene = items[5]; sqanti_assoc_trans = items[6]; strand = items[7]
		#print(newID+'\t'+str(exons))
		if newID == 'PB.7697.1' or newID == 'PB.7374.1_dup2' or newID == 'PB.7374.1':
			continue
		if exons == 1 and coding == 'non_coding':
			continue
		if flair_scaf.startswith('gene'):			
			ncbi_gene_dirn = ncbi_gtf_dic[flair_scaf][0]
			if newID in verified_dic.keys():
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,verified_dic[newID])) #; del sub_sqanti_trans_dic[newID]
			elif ncbi_gene_dirn == strand:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,flair_scaf)) #; del sub_sqanti_trans_dic[newID]
			else:
				print("Flair gene starts with gene but strand is different " + newID)
		
		elif structure == 'full-splice_match' and flair_scaf != sqanti_assoc_gene and not flair_scaf.startswith('gene'):			
			ncbi_gene_dirn = ncbi_gtf_dic[sqanti_assoc_gene][0]
			exons_dic1 = {newID:sqanti_trans_dic[newID]}
			exons_dic2 = {sqanti_assoc_trans:ncbi_trans_dic[sqanti_assoc_trans]}
			#if return_percentage_overlap( sqanti_trans_dic[newID],ncbi_trans_dic[sqanti_assoc_trans] ) > 25 and ncbi_gene_dirn == strand:
			if return_percentage_overlap( exons_dic1,exons_dic2 ) > 25 and ncbi_gene_dirn == strand:
				for each_item in sqanti_lines[newID]:
					newfile.write(each_item)
			else:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID))

			#if ncbi_gene_dirn == strand:
			#	if return_percentage_overlap( sqanti_trans_dic[newID],ncbi_trans_dic[sqanti_assoc_trans] ) > 25:
			#		for each_item in sqanti_lines[newID]:
			#			newfile.write(each_item)
			#	elif return_percentage_overlap( sqanti_trans_dic[newID],ncbi_trans_dic[sqanti_assoc_trans] ) < 25:
			#		for each_item in sqanti_lines[newID]:
			#			sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
			#			newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID))						
			#if sqanti_assoc_trans not in ncbi_CDS_dic.keys(): #Check to see if associated transcript is not coding
			#	if ncbi_gene_dirn == strand:
			#		if return_percentage_overlap() < 25:
			#			
			#		for each_item in sqanti_lines[newID]:
			#			newfile.write(each_item); del sub_sqanti_trans_dic[newID] #If the associated gene does not have a CDS just add gene ID
			#	else:
			#		print("FSM but strand is different " + newID)
			#elif sqanti_assoc_trans in ncbi_CDS_dic.keys() and ncbi_gene_dirn == strand: #
			#	overlap = find_complete_overlapping_exons(strand,sqanti_trans_dic[newID],ncbi_CDS_dic[sqanti_assoc_trans])
			#	if len(overlap) > 0:
			#		for each_item in sqanti_lines[newID]:
			#			newfile.write(each_item); del sub_sqanti_trans_dic[newID] #If associated gene overlaps a CDS we we give it that gene ID
			#	else:
			#		sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
			#		for each_item in sqanti_lines[newID]:
			#			newfile.write(each_item.replace(sqanti_quasi_gene_ID,tofu_gene_ID)); del sub_sqanti_trans_dic[newID]
			#elif sqanti_assoc_trans in ncbi_CDS_dic.keys() and ncbi_gene_dirn != strand:
			#	sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
			#	for each_item in sqanti_lines[newID]:
			#		newfile.write(each_item.replace(sqanti_quasi_gene_ID,tofu_gene_ID)); del sub_sqanti_trans_dic[newID]

		elif structure == 'incomplete-splice_match' and flair_scaf != sqanti_assoc_gene and not flair_scaf.startswith('gene'):			
			#seems quite easy, if multi-exon give it sqanti associated gene ID, if not just assign tofu gene ID. Forget overlaps coz we cannot call them genes and direction does not matter
			ncbi_gene_dirn = ncbi_gtf_dic[sqanti_assoc_gene][0]
			if exons > 1 and ncbi_gene_dirn == strand:
				for each_item in sqanti_lines[newID]:
					newfile.write(each_item) #; del sub_sqanti_trans_dic[newID]
			else:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID)) #; del sub_sqanti_trans_dic[newID]
		
		elif structure == 'novel_not_in_catalog' and flair_scaf != sqanti_assoc_gene and not flair_scaf.startswith('gene'):			
			ncbi_gene_dirn = ncbi_gtf_dic[sqanti_assoc_gene][0]
			if ncbi_gene_dirn == strand:
				for each_item in sqanti_lines[newID]:
					newfile.write(each_item)
			else:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID))
			
			#if sqanti_assoc_trans in ncbi_CDS_dic.keys() or ncbi_gene_dirn != strand:
			#	sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
			#	for each_item in sqanti_lines[newID]:
			#		newfile.write(each_item.replace(sqanti_quasi_gene_ID,tofu_gene_ID)); del sub_sqanti_trans_dic[newID]
			#elif sqanti_assoc_trans not in ncbi_CDS_dic.keys() and ncbi_gene_dirn == strand:
			#	sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
			#	for each_item in sqanti_lines[newID]:
			#		newfile.write(each_item.replace(sqanti_quasi_gene_ID,tofu_gene_ID)); del sub_sqanti_trans_dic[newID]
		
		elif structure == 'novel_in_catalog' and flair_scaf != sqanti_assoc_gene and not flair_scaf.startswith('gene'):
			#If multi exonic and same direction give sqanti associated gene, if single exon but covers a huge chuck of the associated gene give that gene ID else give own tofu transcript ID as gene ID
			ncbi_gene_dirn = ncbi_gtf_dic[sqanti_assoc_gene][0]
			if exons > 1 and ncbi_gene_dirn == strand:
				for each_item in sqanti_lines[newID]:
					newfile.write(each_item)
			elif exons == 1 and ncbi_gene_dirn == strand:
				psuedo_sqanti_assoc_trans = ncbi_gene_trans_dic[sqanti_assoc_gene][0] #get all transcripts for the associated gene and just take the first one in the dictionary, which is randomly picked
				exons_dic1 = {newID:sqanti_trans_dic[newID]}
				exons_dic2 = {psuedo_sqanti_assoc_trans:ncbi_trans_dic[psuedo_sqanti_assoc_trans]}
				if return_percentage_overlap2( exons_dic1,exons_dic2 ) >= 25:
					for each_item in sqanti_lines[newID]:
						newfile.write(each_item)
						#print(psuedo_sqanti_assoc_trans+'\t'+newID + '\t' + str(return_percentage_overlap2( exons_dic1,exons_dic2 ))+'\tOVER')						
				else:
					for each_item in sqanti_lines[newID]:
						sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
						newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID))
						#print(psuedo_sqanti_assoc_trans+'\t'+newID + '\t' + str(return_percentage_overlap2( exons_dic1,exons_dic2 ))+'\tUNDER')
			else:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID))
		
		elif structure == 'genic':
			ncbi_gene_dirn = ncbi_gtf_dic[sqanti_assoc_gene.split('_')[0]][0]
			psuedo_sqanti_assoc_trans = ncbi_gene_trans_dic[sqanti_assoc_gene.split('_')[0]][0]
			exons_dic1 = {newID:sqanti_trans_dic[newID]}
			exons_dic2 = {psuedo_sqanti_assoc_trans:ncbi_trans_dic[psuedo_sqanti_assoc_trans]}
			if newID in genic_verified.keys():
				psuedo_sqanti_assoc_trans = ncbi_gene_trans_dic[genic_verified[newID]][0]
			overlap_perc_x = return_percentage_overlap2( exons_dic1,exons_dic2 )
			if overlap_perc_x >= 25 and ncbi_gene_dirn == strand:
				for each_item in sqanti_lines[newID]:
					newfile.write(each_item)
			elif overlap_perc_x >= 15 and ncbi_gene_dirn == strand and 'CDS' in sqanti_transCDS_dic[newID]:
				for each_item in sqanti_lines[newID]:
					newfile.write(each_item)
			else:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID))
			#print(psuedo_sqanti_assoc_trans+'\t'+newID + '\t' + str(return_percentage_overlap2( exons_dic1,exons_dic2 ))+'\tOVER')
			
		elif structure == 'fusion' and not flair_scaf.startswith('gene'):
			if newID in verified_dic.keys(): 
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,verified_dic[newID]))
			else:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,tofu_gene_ID))
					
		elif structure == 'intergenic' or structure == 'antisense':
			if newID == 'PB.9001.17_dup2' or newID == 'PB.7374.1_dup2':
				sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
				newfile.write(each_item.replace(sqanti_quasi_gene_ID,tofu_gene_IDx))
				continue
			#print("Processing " + newID) #; print(sqanti_gene_dic[tofu_gene_IDx])
			#find overlapping intergenic and assign tofu gene I
			if len(sqanti_gene_dic[tofu_gene_IDx]) == 1:
				for each_item in sqanti_lines[newID]:
					if not gene_IDPat.search(each_item):
						newfile.write(each_item.replace('""', '"%s"' %tofu_gene_IDx))
					else:
						sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
						newfile.write(each_item.replace(sqanti_quasi_gene_ID,tofu_gene_IDx))
			elif len(sqanti_gene_dic[tofu_gene_IDx]) > 1:
				all_trans_exon_dic = {}
				for trans_i in sqanti_gene_dic[tofu_gene_IDx]:
					if structure == classn_dic[trans_i][4]:
						all_trans_exon_dic[trans_i] = sqanti_trans_dic[trans_i]
				#print(all_trans_exon_dic)
				if len(all_trans_exon_dic.keys()) == 1 and newID in all_trans_exon_dic.keys(): #Sometimes there's this transcript among many that is the only allocated to a certain category. It means we just print it with its own gene ID
					for each_item in sqanti_lines[newID]:
						sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
						newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID))
				else:
					for trans_x,perc in find_overlapping_trans(all_trans_exon_dic).items():
						trans_x1 = trans_x.split('_')[0]
						trans_x2 = trans_x.split('_')[1]
						if perc >= 50:
							combine = sqanti_lines[trans_x1] + sqanti_lines[trans_x2]
							for each_item in combine:
								gene_IDPat = re.compile(r'(gene_id ".*"); t'); sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
								trans_IDx = trans_IDPat.search(each_item).group(1)
								if trans_IDx not in significant_overs.keys():
									significant_overs.setdefault(trans_IDx,[each_item.replace(sqanti_quasi_gene_ID,'gene_id "%s"' %tofu_gene_IDx)])
								elif trans_IDx in significant_overs.keys():
									significant_overs[trans_IDx].append(each_item.replace(sqanti_quasi_gene_ID,'gene_id "%s"' %tofu_gene_IDx))
							gene_IDPat = re.compile(r'gene_id "([^"]+)"')
						else:
							combine = sqanti_lines[trans_x1] + sqanti_lines[trans_x2]
							#print('Combine ...'); print(combine)
							for each_item in combine:
								gene_IDPat = re.compile(r'(gene_id ".*"); t'); sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
								trans_IDx = trans_IDPat.search(each_item).group(1)
								if trans_IDx not in all_overs.keys():
									all_overs.setdefault(trans_IDx,[each_item.replace(sqanti_quasi_gene_ID,'gene_id "%s"' %trans_IDx)])
								elif trans_IDx in all_overs.keys():
									all_overs[trans_IDx].append(each_item.replace(sqanti_quasi_gene_ID,'gene_id "%s"' %trans_IDx))
							gene_IDPat = re.compile(r'gene_id "([^"]+)"')
	#print(all_overs)
	for id_1,listx in significant_overs.items():
		for each_item in set(listx):
			newfile.write(each_item)
	for id_1,listx in all_overs.items():
		if id_1 not in significant_overs.keys():
			for each_item in set(listx):
				newfile.write(each_item)
	newfile.close()

def sqanti_gff_edit2(sqanti_gtf, ncbi_gtf, classn):
	dirc = os.path.dirname(sqanti_gtf)
	newfile = open(sqanti_gtf+'_updated', 'w')
	gene_IDPat = re.compile(r'gene_id "([^"]+)"')
	trans_IDPat = re.compile(r'transcript_id "([^"]+)"') 
	ncbi_gtf_dic = {}
	ncbi_trans_dic = {}
	ncbi_CDS_dic = {}
	ncbi_gene_trans_dic = {}
	a=b=c=d=e=f=g=h=i =0
	for line in open(ncbi_gtf, 'r'):
		if len(line.split()) == 10:
			gene_ID = trans_ID = trans_IDPat.search(line).group(1)
		else:
			gene_ID = gene_IDPat.search(line).group(1)	   
			trans_ID = trans_IDPat.search(line).group(1)
		exon_start = min( line.split()[3],line.split()[4] )
		exon_end = max( line.split()[3],line.split()[4] )
		if gene_ID not in ncbi_gene_trans_dic.keys() and line.split()[2] == 'exon':
			ncbi_gene_trans_dic.setdefault(gene_ID,[trans_ID])
		elif gene_ID in ncbi_gene_trans_dic.keys() and line.split()[2] == 'exon' and trans_ID not in ncbi_gene_trans_dic[gene_ID]:
			ncbi_gene_trans_dic[gene_ID].append(trans_ID)
		if gene_ID not in ncbi_gtf_dic.keys():
			ncbi_gtf_dic.setdefault(gene_ID,[line.split()[6],line.split()[2]])
		elif gene_ID in ncbi_gtf_dic.keys() and line.split()[2] not in ncbi_gtf_dic[gene_ID]:
			ncbi_gtf_dic[gene_ID].append(line.split()[2])
		if trans_ID not in ncbi_trans_dic.keys() and line.split()[2] == 'exon':
			ncbi_trans_dic.setdefault(trans_ID,[exon_start+'_'+exon_end])
		elif trans_ID in ncbi_trans_dic.keys() and line.split()[2] == 'exon':
			ncbi_trans_dic[trans_ID].append(exon_start+'_'+exon_end)
		elif trans_ID not in ncbi_CDS_dic.keys() and line.split()[2] == 'CDS':
			ncbi_CDS_dic.setdefault(trans_ID,[exon_start+'_'+exon_end])
		elif trans_ID in ncbi_CDS_dic.keys() and line.split()[2] == 'CDS':
			ncbi_CDS_dic[trans_ID].append(exon_start+'_'+exon_end)
	
	classn_dic = {}
	for line in open(classn, 'r'):
		items = line.split()
		if items[0] not in classn_dic.keys() and a > 0:
			newID = items[0]; oldID = items[0]; flair_scaf = 'YYY'; coding = items[29]; exons = items[4]; structure = items[5]; sqanti_assoc_gene = items[6]; sqanti_assoc_trans = items[7]; strand = items[2]
			classn_dic.setdefault(newID,[ oldID, flair_scaf, coding, exons, structure, sqanti_assoc_gene, sqanti_assoc_trans, strand ])
		a += 1
	
	sqanti_gene_dic = {}
	sqanti_trans_dic = {}
	sqanti_transCDS_dic = {}
	sqanti_lines = {}
	for line in open(sqanti_gtf,'r'):
		gene_IDPat = re.compile(r'gene_id "([^"]+)"')
		trans_ID = trans_IDPat.search(line).group(1)
		tofu_gene_ID = '.'.join(trans_ID.split('.')[:2]) #trans_ID.replace('.'+trans_ID.split('.')[-1],'')
		exon_start = min( line.split()[3],line.split()[4] )
		exon_end = max( line.split()[3],line.split()[4] )
		#add all gff lines to dic
		if trans_ID not in sqanti_lines.keys():
			sqanti_lines.setdefault(trans_ID,[line])
		elif trans_ID in sqanti_lines.keys():
			sqanti_lines[trans_ID].append(line)
		#add all tofu genes and their associated transcripts according to tofu
		if tofu_gene_ID not in sqanti_gene_dic.keys():
			sqanti_gene_dic.setdefault(tofu_gene_ID,[trans_ID])
		elif tofu_gene_ID in sqanti_gene_dic.keys() and trans_ID not in sqanti_gene_dic[tofu_gene_ID]:
			sqanti_gene_dic[tofu_gene_ID].append(trans_ID)
		#add all tofu transcripts and their associated exon coordinates to dic
		if trans_ID not in sqanti_trans_dic.keys() and line.split()[2] == 'exon':
			sqanti_trans_dic.setdefault(trans_ID,[exon_start+'_'+exon_end])
		elif trans_ID in sqanti_trans_dic.keys() and line.split()[2] == 'exon':
			sqanti_trans_dic[trans_ID].append(exon_start+'_'+exon_end)
		#add tofu transcripts and their features whether transcript,exon,CDS
		if trans_ID not in sqanti_transCDS_dic.keys():
			sqanti_transCDS_dic.setdefault(trans_ID,[line.split()[2]])
		elif trans_ID in sqanti_transCDS_dic.keys() and line.split()[2] not in sqanti_transCDS_dic[trans_ID]:
			sqanti_transCDS_dic[trans_ID].append(line.split()[2])
	sub_sqanti_trans_dic = sqanti_trans_dic	
	
	all_overs = {} #; print(sqanti_gene_dic); print('No no no')
	significant_overs = {}; genic_verified = {}; verified_dic = {}
	
	for newID,items in classn_dic.items():
		tofu_gene_IDx = '.'.join(newID.split('.')[:2]) #newID.replace('.'+newID.split('.')[-1],'')
		oldID = items[0]; flair_scaf = items[1]; coding = items[2]; exons = int(items[3]); structure = items[4]; sqanti_assoc_gene = items[5]; sqanti_assoc_trans = items[6]; strand = items[7]
		newID_ed = newID.replace('.'+newID.split('.')[-1], '') #'.'.join(newID.split('.')[0:2]) #Here, we need to remove the last transcript identifier so that we have only gene ID left.
		#print(newID+'\t'+str(exons))
		if exons == 1 and coding == 'non_coding':
			continue
		
		elif structure == 'full-splice_match':			
			ncbi_gene_dirn = ncbi_gtf_dic[sqanti_assoc_gene][0]
			exons_dic1 = {newID:sqanti_trans_dic[newID]}
			exons_dic2 = {sqanti_assoc_trans:ncbi_trans_dic[sqanti_assoc_trans]}
			#if return_percentage_overlap( sqanti_trans_dic[newID],ncbi_trans_dic[sqanti_assoc_trans] ) > 25 and ncbi_gene_dirn == strand:
			if return_percentage_overlap( exons_dic1,exons_dic2 ) >= 25 and ncbi_gene_dirn == strand:
				for each_item in sqanti_lines[newID]:
					newfile.write(each_item)
			else:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID_ed)) 

		elif structure == 'incomplete-splice_match':			
			#seems quite easy, if multi-exon give it sqanti associated gene ID, if not just assign tofu gene ID. Forget overlaps coz we cannot call them genes and direction does not matter
			ncbi_gene_dirn = ncbi_gtf_dic[sqanti_assoc_gene][0]
			if exons > 1 and ncbi_gene_dirn == strand:
				for each_item in sqanti_lines[newID]:
					newfile.write(each_item) #; del sub_sqanti_trans_dic[newID]
			else:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID_ed)) #; del sub_sqanti_trans_dic[newID]
		
		elif structure == 'novel_not_in_catalog' :			
			ncbi_gene_dirn = ncbi_gtf_dic[sqanti_assoc_gene][0]
			if ncbi_gene_dirn == strand:
				for each_item in sqanti_lines[newID]:
					newfile.write(each_item)
			else:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID_ed))
		
		elif structure == 'novel_in_catalog':
			#If multi exonic and same direction give sqanti associated gene, if single exon but covers a huge chuck of the associated gene give that gene ID else give own tofu transcript ID as gene ID
			ncbi_gene_dirn = ncbi_gtf_dic[sqanti_assoc_gene][0]
			if exons > 1 and ncbi_gene_dirn == strand:
				for each_item in sqanti_lines[newID]:
					newfile.write(each_item)
			elif exons == 1 and ncbi_gene_dirn == strand:
				psuedo_sqanti_assoc_trans = ncbi_gene_trans_dic[sqanti_assoc_gene][0] #get all transcripts for the associated gene and just take the first one in the dictionary, which is randomly picked
				exons_dic1 = {newID:sqanti_trans_dic[newID]}
				exons_dic2 = {psuedo_sqanti_assoc_trans:ncbi_trans_dic[psuedo_sqanti_assoc_trans]}
				if return_percentage_overlap2( exons_dic1,exons_dic2 ) >= 25:
					for each_item in sqanti_lines[newID]:
						newfile.write(each_item)
						#print(psuedo_sqanti_assoc_trans+'\t'+newID + '\t' + str(return_percentage_overlap2( exons_dic1,exons_dic2 ))+'\tOVER')						
				else:
					for each_item in sqanti_lines[newID]:
						sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
						newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID_ed))
						#print(psuedo_sqanti_assoc_trans+'\t'+newID + '\t' + str(return_percentage_overlap2( exons_dic1,exons_dic2 ))+'\tUNDER')
			else:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID_ed))
		
		elif structure == 'genic':
			ncbi_gene_dirn = ncbi_gtf_dic[sqanti_assoc_gene.split('_')[0]][0]
			psuedo_sqanti_assoc_trans = ncbi_gene_trans_dic[sqanti_assoc_gene.split('_')[0]][0]
			exons_dic1 = {newID:sqanti_trans_dic[newID]}
			exons_dic2 = {psuedo_sqanti_assoc_trans:ncbi_trans_dic[psuedo_sqanti_assoc_trans]}
			if newID in genic_verified.keys():
				psuedo_sqanti_assoc_trans = ncbi_gene_trans_dic[genic_verified[newID]][0]
			overlap_perc_x = return_percentage_overlap2( exons_dic1,exons_dic2 )
			if overlap_perc_x >= 25 and ncbi_gene_dirn == strand:
				for each_item in sqanti_lines[newID]:
					newfile.write(each_item)
			elif overlap_perc_x >= 15 and ncbi_gene_dirn == strand and 'CDS' in sqanti_transCDS_dic[newID]:
				for each_item in sqanti_lines[newID]:
					newfile.write(each_item)
			else:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID_ed))
			#print(psuedo_sqanti_assoc_trans+'\t'+newID + '\t' + str(return_percentage_overlap2( exons_dic1,exons_dic2 ))+'\tOVER')
			
		elif structure == 'fusion':
			if newID in verified_dic.keys(): 
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,verified_dic[newID]))
			else:
				for each_item in sqanti_lines[newID]:
					sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
					newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID_ed))
					
		elif structure == 'intergenic' or structure == 'antisense':
			if newID == 'PB.9001.17_dup2' or newID == 'PB.7374.1_dup2':
				sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
				newfile.write(each_item.replace(sqanti_quasi_gene_ID,tofu_gene_IDx))
				continue
			#print("Processing " + newID) #; print(sqanti_gene_dic[tofu_gene_IDx])
			#find overlapping intergenic and assign tofu gene I
			if len(sqanti_gene_dic[tofu_gene_IDx]) == 1:
				for each_item in sqanti_lines[newID]:
					if not gene_IDPat.search(each_item):
						newfile.write(each_item.replace('""', '"%s"' %tofu_gene_IDx))
					else:
						sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
						newfile.write(each_item.replace(sqanti_quasi_gene_ID,tofu_gene_IDx))
			elif len(sqanti_gene_dic[tofu_gene_IDx]) > 1:
				all_trans_exon_dic = {}
				for trans_i in sqanti_gene_dic[tofu_gene_IDx]:
					if structure == classn_dic[trans_i][4]:
						all_trans_exon_dic[trans_i] = sqanti_trans_dic[trans_i]
				#print(all_trans_exon_dic)
				if len(all_trans_exon_dic.keys()) == 1 and newID in all_trans_exon_dic.keys(): #Sometimes there's this transcript among many that is the only allocated to a certain category. It means we just print it with its own gene ID
					for each_item in sqanti_lines[newID]:
						sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
						newfile.write(each_item.replace(sqanti_quasi_gene_ID,newID))
				else:
					for trans_x,perc in find_overlapping_trans(all_trans_exon_dic).items():
						trans_x1 = trans_x.split('_')[0]
						trans_x2 = trans_x.split('_')[1]
						if perc >= 50:
							combine = sqanti_lines[trans_x1] + sqanti_lines[trans_x2]
							for each_item in combine:
								gene_IDPat = re.compile(r'(gene_id ".*"); t'); sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
								trans_IDx = trans_IDPat.search(each_item).group(1)
								if trans_IDx not in significant_overs.keys():
									significant_overs.setdefault(trans_IDx,[each_item.replace(sqanti_quasi_gene_ID,'gene_id "%s"' %tofu_gene_IDx)])
								elif trans_IDx in significant_overs.keys():
									significant_overs[trans_IDx].append(each_item.replace(sqanti_quasi_gene_ID,'gene_id "%s"' %tofu_gene_IDx))
							gene_IDPat = re.compile(r'gene_id "([^"]+)"')
						else:
							combine = sqanti_lines[trans_x1] + sqanti_lines[trans_x2]
							#print('Combine ...'); print(combine)
							for each_item in combine:
								gene_IDPat = re.compile(r'(gene_id ".*"); t'); sqanti_quasi_gene_ID = gene_IDPat.search(each_item).group(1)
								trans_IDx = trans_IDPat.search(each_item).group(1)
								if trans_IDx not in all_overs.keys():
									all_overs.setdefault(trans_IDx,[each_item.replace(sqanti_quasi_gene_ID,'gene_id "%s"' %trans_IDx)])
								elif trans_IDx in all_overs.keys():
									all_overs[trans_IDx].append(each_item.replace(sqanti_quasi_gene_ID,'gene_id "%s"' %trans_IDx))
							gene_IDPat = re.compile(r'gene_id "([^"]+)"')
	#print(all_overs)
	for id_1,listx in significant_overs.items():
		for each_item in set(listx):
			newfile.write(each_item)
	for id_1,listx in all_overs.items():
		if id_1 not in significant_overs.keys():
			for each_item in set(listx):
				newfile.write(each_item)
	newfile.close()
	
def find_overlapping_trans(exons_dic): #{PB.12.1:[12_25,30_40,56_70],PB.12.2:[12_25,30_40,56_70]}
	#print('Yes'); print(exons_dic)
	overlapping = {}
	new_dic = {}
	for trans,exons in exons_dic.items():
		new_dic[trans] = exons
	for trans,exons in exons_dic.items(): #We loop over each transcript and its associated exons		
		del new_dic[trans] #we create a new dictionary that eliminates the transcript we are examining
		#print(new_dic)
		if len(new_dic.keys()) > 0:
			for x,y in new_dic.items():
				overlap_perc = return_percentage_overlap2({trans:exons},{x:y})
				overlapping[trans+'_'+x] = overlap_perc
	#print(overlapping)
	return(overlapping)

def find_overlapping_ranges(file1, file2, file3):
	#print('Yes'); print(exons_dic)
	file1_dic = {}
	file2_dic = {}
	verison = ''
	if 'V3' in file1:
		version = '_V3'
	elif 'V4' in file1:
		version = '_V4'
	new_file = open(file3+version+'_overlaps', 'w')
	for line in open(file1, 'r'):
		line_items = line.split()
		file1_dic[line_items[0]] = line_items[8]+':'+line_items[9]
	for line in open(file2, 'r'):
		line_items = line.split()
		file2_dic[line_items[0]] = line_items[8]+':'+line_items[9]
			
	for line in open(file3, 'r'): #We loop over each transcript and its associated exons		
		mutn = line.split()[0]
		for prim,rang in file1_dic.items():
			if int(mutn) in range(int(rang.split(':')[0]),int(rang.split(':')[1])+1):
				new_file.write('\t'.join([prim,rang,mutn+'\n']))
		for prim,rang in file2_dic.items():
			if int(mutn) in range(int(rang.split(':')[0]),int(rang.split(':')[1])+1):
				new_file.write('\t'.join([prim,rang,mutn+'\n']))
	new_file.close()

#def find_complete_overlapping_exons(dirn,exons_dic1,exons_dic2): #{PB.12.1:[12_25,30_40,56_70],PB.12.2:[12_25,30_40,56_70]}
#	overlapping = []
#	for trans,exons in exons_dic1.items(): #We loop over each transcript and its associated exons
#		#new_dic = exons_dic 
#		#del new_dic[trans] #we create a new dictionary that eliminates the transcript we are examining
#		for each_exon in exons: #We loop over all exons for this transcript
#			left = int(each_exon.split('_')[0])
#			right = int(each_exon.split('_')[1])
#			for x,y in exons_dic2: #Now we loop over each of the transcript in the new dictionary
#				for n in y: #We examin their exons
#					nleft = int(n.split('_')[0])
#					nright = int(n.split('_')[1])
#					if nleft in range(left,right) and nright in range(left,right) and dirn == '+':
#						overlapping += [trans+':'+x]
#					elif nleft in range(left,right) and right in range(nleft,nright) and dirn == '+':
#						overlapping += [trans+':'+x]						
#					elif nright in range(left,right) and right in range(nleft,nright) and dirn == '-':
#						overlapping += [trans+':'+x]
#	return(overlapping)
	
def return_percentage_overlap(exons_dic1,exons_dic2): #{PB.12.1:[12_25,30_40,56_70],PB.12.2:[12_25,30_40,56_70]}
	overlapping = 0
	for trans,exons in exons_dic1.items(): #We loop over each transcript and its associated exons
		for each_exon in exons: #We loop over all exons for this transcript
			left = int(each_exon.split('_')[0])
			right = int(each_exon.split('_')[1])
			for x,y in exons_dic2.items(): #Now we loop over each of the ref gene in the new dictionary
				for n in y: #We examin their exons
					nleft = int(n.split('_')[0])
					nright = int(n.split('_')[1])
					denom = nright - nleft
					if left not in range(nleft,nright) and right in range(nleft,nright): # and dirn == '+':
						overlapping = round(100 * (right - nleft)/denom)
					elif left in range(nleft,nright) and right in range(nleft,nright): # and dirn == '+':
						overlapping = round(100 * (right - left)/denom)
					elif left in range(nleft,nright) and right not in range(nleft,nright): # and dirn == '-':
						overlapping = round(100 * (nright - left)/denom)
					elif nleft in range(left,right) and nright in range(left,right): # and dirn == '-':
						overlapping = 100
	return(overlapping)

def return_percentage_overlap2(exons_dic1,exons_dic2): #{PB.12.1:[12_25,30_40,56_70],PB.12.2:[12_25,30_40,56_70]}
	#print('YES')
	overlapping = 0
	for trans,exons in exons_dic1.items(): #We loop over each transcript and its associated exons
		if len(exons) > 1:
			exons = [ '_'.join([str(min(int(exons[0].split('_')[0]),int(exons[0].split('_')[1]))),str(max(int(exons[-1].split('_')[0]),int(exons[-1].split('_')[1])))]) ]
		for each_exon in exons: #We loop over all exons for this transcript
			left = int(each_exon.split('_')[0])
			right = int(each_exon.split('_')[1])
			for x,y in exons_dic2.items(): #Now we loop over each of the ref gene in the new dictionary
				if len(y) > 1:					
					intron_len = x = 0
					while x < (len(y)-1): #for i in range(len(y)):
						intron_len += int(y[x+1].split('_')[0]) - int(y[x].split('_')[1])
						x += 1
					y = [ '_'.join([str(min(int(y[0].split('_')[0]),int(y[0].split('_')[1]))),str(max(int(y[-1].split('_')[0]),int(y[-1].split('_')[1])))]) ]
				else:
					intron_len = 0
				for n in y: #We examin their exons
					nleft = int(n.split('_')[0])
					nright = int(n.split('_')[1])
					denom = nright - nleft - intron_len
					if left not in range(nleft,nright) and right in range(nleft,nright): # and dirn == '+':
						overlapping = round(100 * (right - nleft)/denom)
					elif left in range(nleft,nright) and right in range(nleft,nright): # and dirn == '+':
						overlapping = round(100 * (right - left)/denom)
					elif left in range(nleft,nright) and right not in range(nleft,nright): # and dirn == '-':
						overlapping = round(100 * (nright - left)/denom)
					elif nleft in range(left,right) and nright in range(left,right): # and dirn == '-':
						overlapping = 100
	return(overlapping)

def man_gft_x4(file1, file2, classn, ncbi_gtf, sqanti_gff):
	newfile = open(sqanti_gff+'_edited','w')
	x=a=b=c=d=e=f=g=h=ii=j=k=l=m = 0
	openfile=open(file1,'r')
	lopenfile = openfile.readlines(); l = len(lopenfile); xxx =[]
	
	flair_gene_dic = {}
	for linex in open(file2,'r'):
		flair_gene_dic[linex.split('_')[0].replace('PB.','PB.FLG.')] = linex.split('_')[1].strip()
	
	sqanti_dic = {}
	for linex in open(classn,'r'):
		items = linex.split()
		gene_ID, sqanti_assoc_gene, strand = items[0], items[6], items[2]
		sqanti_dic[gene_ID] = [gene_ID, sqanti_assoc_gene, strand ]
		
	gene_IDPat = re.compile(r'gene_id "([^"]*)"')
	trans_IDPat = re.compile(r'transcript_id "([^"]+)"') 
	ncbi_gtf_dic = {}	
	for line in open(ncbi_gtf, 'r'):
		if len(line.split()) == 10:
			gene_ID = trans_ID = trans_IDPat.search(line).group(1)
		else:
			gene_ID = gene_IDPat.search(line).group(1)	   
			trans_ID = trans_IDPat.search(line).group(1)

		if gene_ID not in ncbi_gtf_dic.keys():
			ncbi_gtf_dic.setdefault(gene_ID,line.split()[6])
			
	verified = {'PB.FLG.96236.2':'gene9302','PB.FLG.29283.3':'PB.FLG.29283.3'}
	length = int(os.popen('wc -l %s' %(sqanti_gff)).read().split()[0])
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(sqanti_gff, 'r'):
		gene_ID = gene_IDPat.search(line).group(1)
		trans_ID = trans_IDPat.search(line).group(1)
		if line.split()[0].startswith('NC_'):
			m+=1
			continue
		if trans_ID+'\n' in lopenfile:
			sqanti_dirn = sqanti_dic[trans_ID][2]
			sqanti_assoc_gene = sqanti_dic[trans_ID][1]
			flair_assoc_gene = flair_gene_dic[trans_ID]
			ncbi_gene_dirn = ncbi_gtf_dic[flair_assoc_gene]
			if sqanti_dirn == ncbi_gene_dirn:
				newfile.write(line.replace(gene_ID,flair_assoc_gene))
				a+=1
			elif sqanti_dirn != ncbi_gene_dirn and len(sqanti_assoc_gene.split('_')) > 1:
				sqanti_associated_genes = sqanti_assoc_gene.split('_')
				sqanti_associated_genes = (x for x in sqanti_associated_genes if x.startswith('gene'))
				b+=1
				if trans_ID in verified.keys():
					newfile.write(line.replace(gene_ID,verified[trans_ID]))
					k +=1
				else:
					for i in sqanti_associated_genes:
						ncbi_gene_dirn = ncbi_gtf_dic[i]
						if sqanti_dirn == ncbi_gene_dirn:
							newfile.write(line.replace(gene_ID,i))
							c+=1
						else:
							newfile.write(line.replace(gene_ID,trans_ID))
							d+=1
			elif sqanti_dirn != ncbi_gene_dirn and len(sqanti_assoc_gene.split('_')) == 1:
				e+=1
				ncbi_gene_dirn = ncbi_gtf_dic[sqanti_assoc_gene]
				if sqanti_dirn == ncbi_gene_dirn:
					newfile.write(line)
					f+=1
				else:
					newfile.write(line.replace(gene_ID,trans_ID))
					g+=1
			else:
				h+=1
			ii += 1
		else:
			j+=1
		x+=1
		bar.update(x)
	bar.finish()
	newfile.close()
	print([a,b,c,d,e,f,g,h,ii,j,x,k,m])
	
def man_gft_x5(gtf, fasta):
	newfile = open(fasta+'_flair_quantify','w')
	x=a=b=c=d=e=f=g=h=ii=j=k=l=m = 0
	
	gene_IDPat = re.compile(r'gene_id "([^"]*)"')
	trans_IDPat = re.compile(r'transcript_id "([^"]+)"') 
	ncbi_gtf_dic = {}	
	for line in open(gtf, 'r'):
		if len(line.split()) == 10:
			gene_ID = trans_ID = trans_IDPat.search(line).group(1)
		else:
			gene_ID = gene_IDPat.search(line).group(1)	   
			trans_ID = trans_IDPat.search(line).group(1)
		if gene_ID.startswith('ERCC'):
			trans_ID = gene_IDPat.search(line).group(1)	   
			gene_ID = trans_IDPat.search(line).group(1)
		if trans_ID not in ncbi_gtf_dic.keys():
			ncbi_gtf_dic.setdefault(trans_ID,trans_ID+'_'+gene_ID+'\n')
		elif trans_ID in ncbi_gtf_dic.keys():
			if trans_ID+'_'+gene_ID+'\n' != ncbi_gtf_dic[trans_ID]:
				print([line,ncbi_gtf_dic[trans_ID]])
				a+=1
	
	#if newID == 'PB.7697.1' or newID == 'PB.101.4':
	#		continue
	lopenfasta = fasta_to_dict(fasta)
	for name,seq in lopenfasta.items():
		name_x = name.strip('>').split()[0].replace('PB.','PB.NEW.')
		if name_x in ncbi_gtf_dic.keys():
			newfile.write('>'+ncbi_gtf_dic[name_x]+seq)
			b += 1
		else:
			c+=1 #; print(name.strip())
	print([a,b,c])
	newfile.close()

def man_gft_x6(gtf):
	newfile1 = open(gtf.replace('_ed.gtf','_ed_Novel.gtf'),'w')
	#newfile2 = open(gtf+'_bad1','w')
	#newfile3 = open(gtf+'_bad4','w')
	#newfile4 = open(gtf+'_bad5','w')
	x=a=b=c=d=e=f=g=h=ii=j=k=l=m = 0

	ref_gene_id_Pat = re.compile(r'ref_gene_id "([^"]*)"')
	reference_id_Pat = re.compile(r'reference_id "([^"]+)"') 
	
	for line in open(gtf, 'r'):
		try:
			ref_gene_id = ref_gene_id_Pat.search(line).group(1)
		except:
			ref_gene_id = ''
		try:
			reference_id = reference_id_Pat.search(line).group(1)
		except:
			reference_id = ''
		
		'''
		if ref_gene_id != '' and reference_id != '':
			continue
		elif ref_gene_id == '' and reference_id == '':
			newfile1.write(line)
		elif ref_gene_id != '' and reference_id == '':
			newfile1.write(line)
			newfile2.write(line)
		elif ref_gene_id == '' and reference_id != '':
			newfile3.write(line)
		else:
			newfile4.write(line)
		'''
		'''
		if ref_gene_id != '' and reference_id != '': #Find all transcripts that StringTie already assigned to a gene and transcript.
			newfile1.write(line)
		'''
		if ref_gene_id == '' and reference_id == '': #Find all transcripts that StringTie has not assigned to a gene and transcript;, these are novel transcripts/genes.
			newfile1.write(line)
	#newfile1.close(); newfile2.close(); newfile3.close(); newfile4.close()
	newfile1.close()
	
def man_gft_x7(gtf):
	newfile1 = open(gtf+'_ed','w')
	newfile2 = open(gtf+'_ed2','w')
	gene_IDPat = re.compile(r'gene_id "([^"]*)"')
	trans_IDPat = re.compile(r'transcript_id "([^"]+)"') 	
	for line in open(gtf, 'r'):
		try:
			gene_ID = gene_IDPat.search(line).group(1)	   
			trans_ID = trans_IDPat.search(line).group(1)
		except:
			newfile2.write(line)
		if gene_ID == trans_ID:
			newfile1.write(line)
	newfile1.close(); newfile2.close()
	
def man_gft_x8(gtf):
	newfile1 = open(gtf+'_ed','w')
	#newfile2 = open(gtf+'_ed2','w')
	gene_IDPat = re.compile(r'gene_id "([^"]*)"')
	trans_IDPat = re.compile(r'transcript_id "([^"]+)"')
	all_transcpts = {}
	transcpt_exon = {}
	scafold = {}
	for line in open(gtf, 'r'):
		scaf = line.split()[0]
		try:
			gene_ID = gene_IDPat.search(line).group(1)	   
			trans_ID = trans_IDPat.search(line).group(1)
		except:
			#newfile2.write(line)
			continue
		if trans_ID not in all_transcpts.keys():
			all_transcpts.setdefault(trans_ID,[line])
		else:
			all_transcpts[trans_ID].append(line)
		if trans_ID not in transcpt_exon.keys():
			if line.split()[2] == 'exon':
				transcpt_exon.setdefault(trans_ID,['exon'])
		else:
			if line.split()[2] == 'exon':
				transcpt_exon[trans_ID].append('exon')
		if scaf not in scafold.keys():
			scafold.setdefault(scaf,[trans_ID])
		else:
			if trans_ID not in scafold[scaf]:
				scafold[scaf].append(trans_ID)
	gd_scafs = []
	gd_trans = []
	for x,y in scafold.items():
		if len(y) > 1:
			gd_scafs += y
	for x,y in transcpt_exon.items():
		if len(y) > 1 and x in gd_scafs and not x.startswith('gene'):
			newfile1.write(''.join(all_transcpts[x])) #gd_trans += all_transcpts[x] #newfile1.write(''.join(all_transcpts[x]))
	newfile1.close()
	'''
	genes = {}
	trans = {}
	for i in gd_trans:
		gene_ID = gene_IDPat.search(line).group(1)	   
		trans_ID = trans_IDPat.search(line).group(1)
		gene_id = trans_ID.replace('.'+trans_ID.split('.')[-1], '')
		if gene_id not in genes.keys():
			genes.setdefault(gene_id,[i])
		else:
			genes[gene_id].append(i)
	for x,y in genes.items():
	'''	
def man_gft_x9(sqanti_gtf, stringtie_gtf, file1, file2):
	newfile1 = open(sqanti_gtf.replace('.gtf','_NOVEL.gtf'),'w')
	newfile2 = open(sqanti_gtf+'_ed2','w')
	gene_IDPat = re.compile(r'gene_id "([^"]*)"')
	trans_IDPat = re.compile(r'transcript_id "([^"]+)"')
	ref_gene_id_Pat = re.compile(r'ref_gene_id "([^"]*)"')
	reference_id_Pat = re.compile(r'reference_id "([^"]+)"')
	nov_trans = []
	known_trans = []
	strg_dic = {}
	for line in open(file1, 'r'):
		nov_trans += [line.strip()]
	for line in open(file2, 'r'):
		try:
			ref_gene_id = ref_gene_id_Pat.search(line).group(1)
		except:
			ref_gene_id = ''
		try:
			reference_id = reference_id_Pat.search(line).group(1)
		except:
			reference_id = ''
		if ref_gene_id != '' and reference_id != '':
			known_trans += [trans_IDPat.search(line).group(1)]
	print(known_trans)
	for line in open(stringtie_gtf, 'r'):
		try:
			gene_ID = gene_IDPat.search(line).group(1)	   
			trans_ID = trans_IDPat.search(line).group(1)
		except:
			newfile2.write(line); gene_ID = ''; trans_ID = ''
		if trans_ID not in strg_dic.keys():
			strg_dic.setdefault(trans_ID,gene_ID)
			
	for line in open(sqanti_gtf, 'r'):
		try:
			gene_ID = gene_IDPat.search(line).group(1)	   
			trans_ID = trans_IDPat.search(line).group(1)
		except:
			newfile2.write(line); gene_ID = ''; trans_ID = ''
		if trans_ID in nov_trans and trans_ID not in known_trans:
			#gene_gene_IDx = '.'.join(trans_ID.split('.')[:2])
			if len(gene_ID.split('_')) > 1 or gene_ID.startswith('novel'):
				newfile1.write(line.replace(gene_ID, strg_dic[trans_ID]))
			else:
				newfile1.write(line)
		else:
			newfile1.write(line)
	newfile1.close(); newfile2.close()
	
def man_gft_x10(gtf):
	newfile1 = open(gtf.replace('.gtf','_ed.gtf','w'))
	newfile2 = open(gtf+'_bad1','w')

	gene_id_Pat = re.compile(r'gene_id "([^"]*)"')
	transcript_id_Pat = re.compile(r'transcript_id "([^"]+)"') 
	
	for line in open(gtf, 'r'):
		try:
			gene_id = gene_id_Pat.search(line).group(1)
		except:
			gene_id = ''
			newfile2.write(line)
		try:
			transcript_id = transcript_id_Pat.search(line).group(1)
		except:
			transcript_id = ''
			newfile2.write(line)
		if gene_id == '':
			newfile1.write(line.strip()+ ' gene_id "%s";\n' %transcript_id)
		else:
			newfile1.write(line)
	newfile1.close(); newfile2.close()
	
def man_gft_x11(gtf1, gtf2, file1):
	newfile1 = open(gtf1.replace('.gtf','_ed.gtf'),'w')
	newfile2 = open(gtf1+'_DELETE','w')
	gene_IDPat = re.compile(r'gene_id "([^"]*)"')
	trans_IDPat = re.compile(r'transcript_id "([^"]+)"')
	filtered_genes = []
	included_genes = []
	a=b=c=d=e= 0
	for line in open(file1, 'r'):
		k = line.strip().split('.')
		m = len(k)
		o = '.'.join( line.strip().split('.')[0:(m-1)] )
		included_genes += [o]
	print(len(included_genes)); print('The total number of included genes is ' + str(len(list(set(included_genes))))); print(included_genes[0:3])
	for line in open(gtf2, 'r'):
		try:
			gene_ID = gene_IDPat.search(line).group(1)	   
			trans_ID = trans_IDPat.search(line).group(1)
		except:
			a += 1; newfile2.write(line)
		filtered_genes += [gene_ID]
	
	new = []
	for i in filtered_genes:
		if i.startswith('ST') or i.startswith('PB'):
			new += [i]
	print('\nThe total number of NOVEL genes is ' + str(len(list(set(new))))); print([x for x in filtered_genes if x.startswith('STRG')][0:3])
	#if 'STRG.30231' in included_genes:
	#	print('Yes')
	#else:
	#	print('No')
	for line in open(gtf1, 'r'):
		try:
			gene_ID = gene_IDPat.search(line).group(1)	   
			trans_ID = trans_IDPat.search(line).group(1)
		except:
			b += 1; newfile2.write(line); continue
		if gene_ID.startswith('ERCC') or gene_ID.startswith('gene') or gene_ID.startswith('rna') or gene_ID in filtered_genes or gene_ID in included_genes:
			newfile1.write(line)
			if gene_ID in filtered_genes and gene_ID not in included_genes:
				c += 1; newfile2.write('*****'+line)			
		else:
			d += 1
	newfile1.close(); newfile2.close()
	print([a,b,c,d])

def man_gft_x12(gtf1):
	gene_IDPat = re.compile(r'gene_id "([^"]*)"')
	new = []
	for line in open(gtf1, 'r'):
		gene_ID = gene_IDPat.search(line).group(1)
		#if gene_ID.startswith('ST') or gene_ID.startswith('PB') or gene_ID.startswith('gene') or gene_ID.startswith('rna') :
		if gene_ID.startswith('ST') or gene_ID.startswith('PB'):
			if gene_ID not in new:
				new += [gene_ID]
	print(len(new))	

def man_gft_x13(gtf1):
	newfile1 = open(gtf1.replace('.gtf','.NovelGenes.gtf'),'w')
	gene_IDPat = re.compile(r'gene_id "([^"]*)"')
	new = []
	for line in open(gtf1, 'r'):
		gene_ID = gene_IDPat.search(line).group(1)
		if gene_ID.startswith('ST') or gene_ID.startswith('PB'):
			if line.split('\t')[2] != "CDS":
				newfile1.write(line)
			if gene_ID not in new:	
				new += [gene_ID]
	newfile1.close(); print(len(new))
	
def man_gft_x14(gtf1):
	dirc = os.path.dirname(gtf1)
	gene_IDPat = re.compile(r'gene_id "([^"]*)"')
	new = []
	for line in open(gtf1, 'r'):
		gene_ID = gene_IDPat.search(line).group(1)
		if gene_ID != '' and gene_ID not in new:
			new += [gene_ID]	
	newfile1 = open(dirc+'/ALL_GENE_NAMES','w')
	newfile1.write('\n'.join(new))
	newfile1.close()
	print(len(new))
	
def average_items(table): #This is a great simple script. It creates a dictionary inside another dictionary and keeps updating that inside dictionary
	isform_dic = {}
	k = open(table, 'r')
	sample_list = k.readlines()[0].strip().split()
	#print(sample_list[0])
	a=0
	#for line in next(open(table, 'r')):
	for line in open(table, 'r'):
		if a == 0:
			a+=1
			continue
		items = line.strip().split()
		isoform = items[0]
		isform_dic.setdefault(isoform,{})
		
		#print(items)
		for i in range(len(items)-1):
			#print(i)
			if sample_list[i+1] not in isform_dic[isoform].keys():
				isform_dic[isoform][sample_list[i+1]] = [items[i+1]]
			else:
				isform_dic[isoform][sample_list[i+1]].append(items[i+1])
		a+=1
	newfile = open(table+'.averages','w')
	
	import pandas as pd
	my_df = pd.DataFrame(sample_list, columns=['sample_list'])
	my_df = list(my_df.drop_duplicates(keep='first')['sample_list'])
	newfile.write('\t'.join(my_df)+'\n')
	for isoformx,samples in isform_dic.items():
		averages = []
		sample_ids2 = []
		for x,y in samples.items():
			averages += [str(eval('+'.join(y))/len(y))]
			sample_ids2 += [x]
		newfile.write(isoformx+'\t'+'\t'.join(averages)+'\n')
	newfile.close()
	
def delete_zxxz(gff3, genes):
	dirc = os.path.dirname(gff3)
	nfile1 = open(dirc+'/'+'reviewer_file', 'w')
	c=d=e=f=g = 0
	items = {}; gene_lines = open(genes, 'r').readlines()
	for line in open(gff3, 'r'):
		line_x = line.split()
		if line_x[0] not in items.keys():
			items.setdefault(line_x[0], [line])
		else:
			items[line_x[0]].append(line)
	for keyx, valuex in items.items():
		if keyx+'\n' in gene_lines or keyx.replace('path1','')+'\n' in gene_lines:
			max_lgFC = [float(x.split()[2]) for x in valuex ]
			nfile1.write(valuex[ max_lgFC.index(max(max_lgFC)) ])
	nfile1.close()
	
def walk_dir1_gunzip(rootDir): 
	#os.system('mkdir %s/pass && mkdir %s/fail' %(rootDir,rootDir))
	#dirs = ['barcode01','barcode02','barcode03','barcode04','barcode05','barcode06','barcode07','barcode08']
	dirs = ['barcode01','barcode02','barcode03','barcode04']
	updir = ''
	for dirName, subdirList, fileList in os.walk(rootDir): #, topdown=False):
		if len(subdirList) == 0 and os.path.basename(dirName) == 'unclassified': # and int(os.path.basename(dirName).strip('barcode')) < 21: #and os.path.basename(dirName) in dirs: #and os.path.basename(dirName) != 'minimap2':
			b_dirName = os.path.basename(dirName)
			print([b_dirName, dirName])
			#'''
			fileList2 = [1,2,3] #fileList2 = [x for x in fileList if x.startswith('barcode')]
			fileList = [x for x in fileList if x.endswith('.fastq.gz')] #fileList = [x for x in fileList if x.endswith('.fastq')]
			if len(fileList) == 0:
				continue #k = 2 #
			n = 0
			bar = progressbar.ProgressBar(maxval=len(fileList), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
			bar.start()
			for filex in fileList:
				if len(fileList2) == 0:
					continue
				else:
					k = 2 #print([b_dirName, dirName])
					#os.system('cd %s && samtools view --threads 23 -b -o %s %s' %(str(dirName), filex.replace('sam','bam'), filex))
					#os.system('cd %s && samtools sort --threads 23 -o %s %s' %(str(dirName), filex.replace('sam','sort.bam'), filex.replace('sam','bam')))
					#os.system('cd %s && samtools index %s' %(str(dirName), filex.replace('sam','sort.bam')))
					#os.system('cd %s && rm %s' %(str(dirName), filex.replace('sam','bam')))
					#os.system('cd %s && /home/banthony/software/Porechop/porechop-runner.py --check_reads 100 -i %s -o %s --threads 23' %(str(dirName),filex,filex.replace('fastq','pchop1.fastq')))
					#os.system('cd %s && bamCoverage --outFileFormat bedgraph --normalizeUsing RPKM -p 20 -b %s -o %s' %(str(dirName),filex,filex.replace('bam','bed')))
				#os.system('cd %s && gunzip -c %s >> %s' %(str(dirName), filex, b_dirName+'.fastq'))
				os.system('gunzip -c %s >> %s' %(str(dirName)+"/"+filex, str(dirName) +'.fastq')) #b_dirName+'.fastq'))
				n += 1
				bar.update(n)
			bar.finish()
			#'''
		#if len(subdirList) == 0 and len(fileList) == 0:
		#	os.system("rm -r %s " %dirName )
				
def walk_dir2_porechop(rootDir):
	print(rootDir); reads_manifest = []
	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0 and 'barcodex' not in dirName:
			b_dirName = os.path.basename(dirName.rstrip('/'))
			fileList2 = [x for x in fileList if x.endswith('.full_length.fq')]
			fileList = [x for x in fileList if x.endswith('.fastq')]
			if len(fileList) == 0:
				continue
			for filex in fileList:
				if len(fileList2) == 0:
					prefix = filex.split('.')[0]
					print([b_dirName, dirName, prefix])
					if b_dirName == 'unclassified_q':
						nn = 2 #os.system('cd %s && /home/banthony/software/Porechop/porechop-runner.py --check_reads 50000 -i %s -o %s --threads 23' %(str(dirName),filex,filex.replace('fastq','pchop1.fastq')))
					else:
						os.system('cd %s && /home/banthony/software/Porechop/porechop-runner.py --check_reads 1000 -i %s -o %s --threads 23' %(str(dirName),filex,filex.replace('fastq','porechop.fastq')))
						#Pychopper module load python/3.9.6 module load hmmer/3.2.1
						'''
						prefix = prefix+'.porechop'
						primers = '/home/banthony/scratch/reads/nanopore/capitata/B002_05_6C2/20220517_2010_2-E11-H11_PAI02157_73737252/fastq_pass/fastq/primers.fasta'
						primer_config = '/home/banthony/scratch/reads/nanopore/capitata/B002_05_6C2/20220517_2010_2-E11-H11_PAI02157_73737252/fastq_pass/fastq/primer_config.txt'
						reads = prefix+'.fastq'
						os.system('cd %s && rm %s.unclassified.fq %s.rescued.fq %s.full_length.fq %s.all.fastq' %(dirName, prefix, prefix, prefix, prefix))
						os.system('cd %s && cdna_classifier.py -t 23 -m edlib -b %s -c %s -r report.pdf -u %s.unclassified.fq -w %s.rescued.fq %s %s.full_length.fq' %(dirName, primers, primer_config, prefix, prefix, reads, prefix))
						run_minimap2_liqa(dirName, prefix)
						'''
						#os.system('cd %s && awk \' {if ($1 ~/^@/) {print ($1\"_FL\")} else print }\' %s.full_length.fq > %s.full_length2.fq' %(dirName, prefix, prefix))
						#os.system('cd %s && awk \' {if ($1 ~/^@/) {print ($1\"_UC\")} else print }\' %s.unclassified.fq > %s.unclassified2.fq' %(dirName, prefix, prefix))
						#os.system('cd %s && awk \' {if ($1 ~/^@/) {print ($1\"_Rs\")} else print }\' %s.rescued.fq > %s.rescued2.fq' %(dirName, prefix, prefix))
						#os.system('cd %s && rm %s.unclassified.fq %s.rescued.fq %s.full_length.fq' %(dirName, prefix, prefix, prefix))
						#os.system('cd %s && mv %s.unclassified2.fq %s.unclassified.fq' %(dirName, prefix, prefix))
						#os.system('cd %s && mv %s.rescued2.fq %s.rescued.fq' %(dirName, prefix, prefix))
						#os.system('cd %s && mv %s.full_length2.fq %s.full_length.fq' %(dirName, prefix, prefix))
						#os.system('cd %s && mv report.pdf %s.report.pdf' %(dirName, prefix))
						#os.system('cd %s && mv cdna_classifier_report.tsv %s.cdna_classifier_report.tsv' %(dirName, prefix))
						
						#os.system("cd %s && cat %s.porechop.full_length.fq %s.porechop.rescued.fq %s.porechop.unclassified.fq > %s.all.porechop.pychop.fq" %(dirName, prefix, prefix, prefix, prefix))
						#reads_manifest += ['\t'.join([prefix,'A','Batch1',dirName+'/'+prefix+'.all.porechop.pychop.fq'])]
						#os.system('cd %s && ~/.local/bin/cutadapt -j 23 -a "A{1000}" -o %s.all.porechop.pychop.cutadapt.fq %s.all.porechop.pychop.fq' %(dirName, prefix, prefix))
						#os.system('cd %s && mv %s.all.porechop.pychop.cutadapt.fasta %s.all.porechop.pychop.cutadapt.fq' %(dirName, prefix, prefix))
						#os.system("cat %s.full_length.fq %s.rescued.fq %s.unclassified.fq | awk '{if(NR%4==1) {printf(\">%s\n\",substr($0,2));} else if(NR%4==12) print;}' > %s.all.fasta" %(dirName, prefix, prefix, prefix, prefix))
				else:
					continue
	#newfile = open('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Plate1_reads_manifest.tsv','w'); newfile.write('\n'.join(reads_manifest)); newfile.close()
	
def run_minimap2_liqa(dirName, prefix):
	#minimap2 export PATH=/home/banthony/software/minimap2-2.24_x64-linux:$PATH module load samtools/1.10 
	folder = re.findall(r'C\d+(\w){1}',prefix)[0]
	folder = prefix.split('.')[0].rstrip(folder)+'H' #folder = prefix.split('.')[0].replace(folder,'H')
	destn = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/'+folder+'/'+prefix.split('.')[0]
	os.makedirs(destn,exist_ok=True)
	#liqa module load python/3.9.6 module load r/4.1.2 
	ref = '/home/banthony/projects/rrg-ioannisr/banthony/analysis/nanopore/capitata/ncbi_genome/genome_ercc/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_plus_EGII_novelgenes.fasta.minimap2.idx'
	unclass = dirName+'/'+prefix+'.unclassified.fq'; resc = dirName+'/'+prefix+'.rescued.fq'; full = dirName+'/'+prefix+'.full_length.fq'
	os.system('cd %s && minimap2 --secondary=no -ax splice -t 23 -o %s.sam %s %s %s %s' %(destn, prefix, ref, unclass, resc, full))
	os.system('cd %s && samtools view --threads 23 -b -o %s.bam %s.sam' %(destn, prefix, prefix))
	os.system('cd %s && samtools sort --threads 23 -o %s.sort.bam %s.bam' %(destn, prefix, prefix))
	os.system('cd %s && samtools index %s.sort.bam' %(destn, prefix))
	os.system('cd %s && rm %s.bam' %(destn, prefix))
	
	#liqa
	gtf = '/home/banthony/projects/rrg-ioannisr/banthony/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_FINAL_flair.gtf'
	os.system('cd %s && liqa -task quantify -refgene %s.refgene -bam %s.sort.bam -out %s -max_distance 10 -f_weight 1' %(destn, gtf, prefix, prefix+'.liqa'))
			  
def walk_dir3(rootDir, details):
	det_dic = {}
	for line in open(details, 'r'):
		barcode = line.strip().split('\t')[-1].replace('BC','')
		subfolder = line.strip().split('\t')[0].split('_')[-2]
		if 'v3' in line:
			folder = 'artic_v3'
		elif 'v4' in line:
			folder = 'artic_v4'
		elif 'Midnight' in line:
			folder = 'midnight'
		elif 'SNAP' in line:
			folder = 'snap'
		det_dic.setdefault(barcode,[folder,subfolder])
		
	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0:
			b_dirName = os.path.basename(dirName); print(b_dirName)
			fileList = [x for x in fileList if x.endswith('.bed')]
			for filex in fileList:
				coverage = 0
				for linex in open(str(dirName)+'/'+filex, 'r'):
					reads = float(linex.strip().split()[-1])
					if float("{:f}".format(reads)) > 10000:
						coverage += 1
				coverage = (100*coverage/600)
				det_dic[filex.split('.')[0].strip('barcode')].append(str(coverage)+'\n')
	newfile = open(rootDir+'/all_coverage.bed', 'w')
	for barcode,details in det_dic.items():
		newfile.write('\t'.join(details))

def read_details(details):
	det_dic = {}
	for line in open(details, 'r'):
		barcode = line.strip().split('\t')[0].strip('BC') #,'')
		subfolder1,subfolder2,subfolder3 = ['','','']
		if len(line.split()) == 2:
			subfolder1 = line.strip().split('\t')[1]
		elif len(line.split()) == 3:
			subfolder1 = line.strip().split('\t')[1].split('_')[-1]; subfolder2 = line.strip().split('\t')[2] #.split('_')[-2]
		elif len(line.split()) == 4:
			subfolder1 = line.strip().split('\t')[1] #.split('_')[-1]
			subfolder2 = line.strip().split('\t')[2] #.split('_')[-2]
			subfolder3 = line.strip().split('\t')[3] #.split('_')[-3]
		if 'v3' in line:
			folder=protocol = 'artic_v3'
		elif 'v4.1' in line:
			folder=protocol = 'artic_v4.1'
		elif 'v4' in line:
			folder=protocol = 'artic_v4'
		elif 'Midnight'.lower() in line:
			folder=protocol = 'midnight'
		elif 'SNAP'.lower() in line:
			folder=protocol = 'snap'
		elif 'Entebbe'.lower() in line:
			folder=protocol = 'entebbe'
		elif 'EBB_' in line:
			folder = 'ebb'; protocol = 'entebbe'
		elif 'Qiaseq' in line:
			folder = protocol = 'Qiaseq'
		else:
			folder=protocol = ''
		#det_dic.setdefault(barcode,[folder,subfolder1,subfolder2])
		det_dic.setdefault(barcode,(protocol+'\t'.join([subfolder1,subfolder2,subfolder3]).replace('\t\t','')).split() ) #det_dic.setdefault(barcode,[protocol,subfolder2,subfolder1])
	return det_dic

def read_details2(details):
	det_dic = {}
	for line in open(details, 'r'):
		barcode,sample_id,sex,plate = line.split()
		barcode = barcode.strip('BC')
		if plate not in det_dic.keys():
			det_dic[plate] = {barcode:[sample_id,sex]}
		else:
			det_dic[plate][barcode] = [sample_id,sex]
	return det_dic

def walk_dir3_count_reads(rootDir, details):
	det_dic = read_details(details)

	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0 and os.path.basename(dirName) != 'unclassified_q':
			b_dirName = os.path.basename(dirName); print([dirName,b_dirName])
			fileList = [x for x in fileList if x.endswith('pchop1.fastq')] #fileList = [x for x in fileList if x.endswith('.pchop1.fastq')]
			for filex in fileList:
				coverage = 0
				length = int(os.popen('wc -l %s' %(str(dirName)+'/'+filex)).read().split()[0])
				reads = (length)/4
				det_dic[filex.split('.')[0].replace('barcode','')].append(str(reads)+'\n')
		elif len(subdirList) == 0 and os.path.basename(dirName) == 'unclassified_q':
			print([dirName,b_dirName])
			fileList = [x for x in fileList if x.endswith('pchop1.fastq')]
			for filex in fileList:
				length = int(os.popen('wc -l %s' %(str(dirName)+'/'+filex)).read().split()[0])
				reads = (length)/4
				det_dic.setdefault('unclassified',['unclassified','.','.',str(reads)+'\n'])
			
	newfile = open(rootDir+'/number_of_pchop_reads.tsv', 'w')
	for barcode,details in det_dic.items():
		newfile.write('\t'.join(details))

def walk_dir3_flagstat(rootDir, details):
	det_dic = read_details(details)
	
	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0 and 'minimap2' not in str(dirName):
			b_dirName = os.path.basename(dirName); print(dirName)
			fileList = [x for x in fileList if x.endswith('.sort.bam')]
			for filex in fileList:
				barcode = filex.split('.')[0].strip('barcode')
				mapped = os.popen("samtools flagstat %s | awk 'NR == 1 || NR == 5{print}' | awk ' {getline a; print($1,a)}' " %(str(dirName)+'/'+filex)).read()
				total = mapped.split()[0]; mapping = mapped.split()[1]; rate = mapped.split()[-3].strip("(").strip('%')
				det_dic[barcode].append('\t'.join([total,mapping,rate,barcode])+'\n')

	newfile = open(rootDir+'/flagstats.tsv', 'w')
	newfile.write('\t'.join(['Protocol','concn','total_reads','mapped_reads','mapping_rate','barcode'])+'\n')
	for barcodex,details in det_dic.items():
		if len('\t'.join(details).split()) > 5:
			newfile.write('\t'.join(details))
			
def walk_dir3_flagstat2(rootDir, details):
	det_dic = {}
	for i in open(details, 'r'):
		det_dic.setdefault(i.split()[1],i.split()[0]+'\t'+i.split()[2])
	newfile = open(rootDir+'/flagstats.tsv', 'w')
	newfile.write('\t'.join(['sample_id2','sample_id','barcode','sex','timepoint','total_reads','mapped_reads','mapping_rate\n']))
	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0 and len([x for x in fileList if x.endswith('.porechop.sam')]) > 0:
			b_dirName = os.path.basename(dirName); print(dirName)
			fileList = [x for x in fileList if x.endswith('.porechop.sam')]
			for filex in fileList:
				sample_id = filex.split('.')[0]
				embryo = re.findall(r'C\d+(\w){1}',sample_id)[0]
				timepoint = sample_id.rstrip(embryo)+'H'
				barcode = det_dic[sample_id].split('\t')[0]; sex = det_dic[sample_id].split('\t')[1]
				sample_id2 = timepoint+barcode+sex
				barcode_details = '\t'.join([sample_id2,sample_id,barcode,sex,timepoint])
				
				mapped = os.popen("samtools flagstat %s | awk 'NR == 1 || NR == 5{print}' | awk ' {getline a; print($1,a)}' " %(str(dirName)+'/'+filex)).read()
				total = mapped.split()[0]; mapping = mapped.split()[1]; rate = mapped.split()[-3].strip("(").strip('%')
				newfile.write('\t'.join([barcode_details,total,mapping,rate])+'\n')
	newfile.close()
	
def walk_dir3_flagstat3(rootDir, details):
	det_dic = {}
	for i in open(details, 'r'):
		det_dic.setdefault(i.split()[1],i.split()[0]+'\t'+i.split()[2])
	#newfile = open('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/transcriptome_flagstats_correct.tsv', 'w')
	newfile = open('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/genome_flagstats_correct.tsv', 'w')
	newfile.write('\t'.join(['sample_id2','sample_id','barcode','sex','timepoint','total_reads','mapped_reads','mapping_rate\n']))
	#newfile.write('\t'.join(['sample_id2','sample_id','barcode','sex',"ERCC","MIT","NOG","GEN","Alignments","errors\n"]))
	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0 and len([x for x in fileList if x.endswith('.all.porechop.pychop.cutadapt.fq')]) > 0:
			b_dirName = os.path.basename(dirName); print(dirName)
			fileList = [x for x in fileList if x.endswith('.all.porechop.pychop.cutadapt.fq')]
			for filex in fileList:
				sample_id = filex.split('.')[0]
				embryo = re.findall(r'C\d+(\w){1}',sample_id)[0]
				timepoint = sample_id.rstrip(embryo)+'H'
				barcode = det_dic[sample_id].split('\t')[0]; sex = det_dic[sample_id].split('\t')[1]
				sample_id2 = timepoint+barcode+sex
				barcode_details = '\t'.join([sample_id2,sample_id,barcode,sex,timepoint])
				
				genome = '/home/banthony/projects/rrg-ioannisr/banthony/analysis/nanopore/capitata/ncbi_genome/genome_ercc/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_plus_EGII_novelgenes.fasta.minimap2.idx'
				#genome = '/home/banthony/projects/rrg-ioannisr/banthony/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/final_files/ERCC_GCF_000347755.3_Ccap_2.1_genomic.transcripts_FINAL_flair_quantify.fasta.minimap2.idx'
				reads = str(dirName)+'/'+filex
				mapped = os.popen(" minimap2 --secondary=no -ax splice -t 23 %s %s | samtools flagstats - | awk 'NR == 1 || NR == 5 {print}' | awk ' {getline a; print($1,a)}' " %(genome, reads)).read()
				total = mapped.split()[0]; mapping = mapped.split()[1]; rate = mapped.split()[-3].strip("(").strip('%')
				results = '\t'.join([total,mapping,rate])
				#mapped = os.popen(" minimap2 -ax splice --secondary=no -t 23 %s %s | awk 'BEGIN{OFS=\"\\t\"} $1 !~/^@/ {if($2 == 4){NOG+=1} else if ($3 ~/^NW_/ || $3 ~/^CAJ/){GEN+=1} else if($3 ~/^NC_/){MIT+=1} else if($3 ~/^ERCC/){ERCC+=1} else {a+=1}; b+=1} END{print(ERCC,MIT,NOG,GEN,b,a)}' " %(genome, reads)).read()
				#ERCC,MIT,NOG,GEN,Alignments,errors = mapped.split('\t'); results = '\t'.join([ERCC,MIT,NOG,GEN,Alignments,errors]).strip()
				newfile.write('\t'.join([barcode_details,results])+'\n')
	newfile.close()
	
def walk_dir3_flagstat4(rootDir, details):
	det_dic = {}
	for i in open(details, 'r'):
		det_dic.setdefault(i.split()[1],i.split()[0]+'\t'+i.split()[2])
	newfile = open(rootDir+'/Cap_all/genome_flagstats.tsv', 'w')
	newfile.write('\t'.join(['sample_id2','sample_id','barcode','sex','timepoint',"ERCC","MIT","NOG","GEN","Alignments",'errors','ERCC_rate','MIT_rate','NOG_rate','GEN_rate','total_reads_genome','mapped_reads_genome','mapping_rate_genome','total_reads_transcriptome','mapped_reads_transcriptome','mapping_rate_transcriptome\n']))
	###
	rd_manfst = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/flair/reads_manifestABC_bulk.txt'
	rd_dic = {}
	for i in open(rd_manfst, 'r'):
		sample_id = i.split('.')[0]
		id1 = sample_id[-1]
		id2 = sample_id.rstrip(id1)
		sample_id = id2.strip('H')+id1
		rd_dic.setdefault(sample_id,i.strip().split('\t')[3])
	###
	os.system("export PATH=/home/banthony/software/minimap2-2.24_x64-linux:$PATH && module load samtools/1.10")
	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0 and len([x for x in fileList if x.endswith('.sort.bam')]) > 0 and 'Cap_all' not in dirName:
			b_dirName = os.path.basename(dirName); print(dirName)
			fileList = [x for x in fileList if x.endswith('.sort.bam')]
			for filex in fileList:
				sample_id = filex.split('.')[0]
				embryo = re.findall(r'C\d+(\w){1}',sample_id)[0]
				timepoint = sample_id.rstrip(embryo)+'H'
				barcode = det_dic[sample_id].split('\t')[0]; sex = det_dic[sample_id].split('\t')[1]
				sample_id2 = timepoint+barcode+sex
				barcode_details = '\t'.join([sample_id2,sample_id,barcode,sex,timepoint])
				#Genome alignment stats
				mapped = os.popen("samtools flagstat %s | awk 'NR == 1 || NR == 5{print}' | awk ' {getline a; print($1,a)}' " %(str(dirName)+'/'+filex)).read()
				total_reads_genome = mapped.split()[0]; mapped_reads_genome = mapped.split()[1]; mapping_rate_genome = mapped.split()[-3].strip("(").strip('%')
				sam_file = str(dirName)+'/'+filex.replace('.sort.bam','.sam')
				mapped = os.popen(" cat %s | awk 'BEGIN{OFS=\"\\t\"} $1 !~/^@/ {if($2 == 4){NOG+=1} else if ($3 ~/^NW_/ || $3 ~/^CAJ/ || $3 ~/^STRG/){GEN+=1} else if($3 ~/^NC_/){MIT+=1} else if($3 ~/^ERCC/){ERCC+=1} else {a+=1}; b+=1} END{print(ERCC,MIT,NOG,GEN,b,a)}' " %(sam_file)).read()
				ERCC,MIT,NOG,GEN,Alignments,errors = mapped.split('\t'); errors = errors if errors.isnumeric() else '0' ; results = '\t'.join([ERCC,MIT,NOG,GEN,Alignments,errors]).strip()
				ERCC_rate,MIT_rate,NOG_rate,GEN_rate = [str(round(100*float(ERCC)/int(Alignments),2)),str(round(100*float(MIT)/int(Alignments),2)),str(round(100*float(NOG)/int(Alignments),2)),str(round(100*float(GEN)/int(Alignments),2))]
				
				#Transcriptome alignment stats
				transcriptome='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_flair_tofu_stringtie.sorted_ed.fa'
				reads=rd_dic[sample_id]
				mapped = os.popen(" minimap2 -ax splice -t 23 %s %s | samtools flagstats - | awk 'NR == 1 || NR == 5 {print}' | awk ' {getline a; print($1,a)}' " %(transcriptome, reads)).read()
				total_reads_transcriptome = mapped.split()[0]; mapped_reads_transcriptome = mapped.split()[1]; mapping_rate_transcriptome = mapped.split()[-3].strip("(").strip('%')
				
				newfile.write('\t'.join([barcode_details,ERCC,MIT,NOG,GEN,Alignments,errors,ERCC_rate,MIT_rate,NOG_rate,GEN_rate,total_reads_genome,mapped_reads_genome,mapping_rate_genome,total_reads_transcriptome,mapped_reads_transcriptome,mapping_rate_transcriptome])+'\n')
				#[ERCC,MIT,NOG,GEN,Alignments,errors,ERCC_rate,MIT_rate,NOG_rate,GEN_rate,total_reads_genome,mapped_reads_genome,mapping_rate_genome,total_reads_transcriptome,mapped_reads_transcriptome,mapping_rate_transcriptome]
	newfile.close()
	
def move_files(directory, details):
	dirc = '/home/banthony/scratch/reads/nanopore/covid/B004/' #os.path.dirname(samfile)
	a=b=c=d=e=f=g=x = 0
	det_dic = {}

	for line in open(details, 'r'):
		barcode = line.strip().split('\t')[-1].replace('BC','')
		subfolder = line.strip().split('\t')[0].split('_')[-2]
		if 'v3' in line:
			folder = 'artic_v3'
		elif 'v4' in line:
			folder = 'artic_v4'
		elif 'Midnight' in line:
			folder = 'midnight'
		elif 'SNAP' in line:
			folder = 'snap'
		det_dic.setdefault(barcode,folder+'/'+subfolder)
		
	files = os.listdir(directory)
	files = [x for x in files if x.startswith('barcode')]
	for filex in files:
		print(filex)
		barcodex = filex.strip('.fastq').replace('barcode','')
		destination = directory + '/' + det_dic[barcodex]
		os.system('cd %s && mv %s %s' %(directory, filex, destination))
		#break
		
def return_folder_x(details):
	det_dic = {}
	'''
	for line in open(details, 'r'):
		barcode = line.strip().split('\t')[1].strip('BC') #.replace('BC','')
		if '.' in line:
			folder,subfolder1 = line.split()[0].split('.')
		else:
			folder = line.split()[0]; subfolder1 = ''
		det_dic.setdefault(barcode,folder+'/'+subfolder1+'/')
	'''
	for line in open(details, 'r'):
		barcode = line.strip().split('\t')[1].strip('BC') #replace('BC','')
		#subfolder1 = line.strip().split('\t')[0].split('_')[-2] #variant
		#subfolder2 = line.strip().split('\t')[0].split('_')[-1].split('.')[0] #level
		if 'v3' in line:
			folder = 'artic_v3'
		elif 'v4.1' in line:
			folder = 'artic_v4.1'
		elif 'Midnight' in line:
			folder = 'midnight'
		elif 'SNAP' in line:
			folder = 'snap'
		elif 'EBB' in line:
			folder = 'ebb'
		elif 'Entebbe' in line:
			folder = 'entebbe'
		elif 'v4L' in line:
			folder = 'artic_v4L'
		elif 'v4C' in line:
			folder = 'artic_v4C'
		else:
			folder = line.strip().split('\t')[0].split('_')[0]
		val = folder #folder+'/'+subfolder1+'/'+subfolder2
		det_dic.setdefault(barcode,val)
	return det_dic
		
def move_files_unzip(rootDir, details):
	a=b=c=d=e=f=g=x = 0
	det_dic = return_folder_x(details)
	
	for dirName, subdirList, fileList in os.walk(rootDir): #, topdown=False):
		base_dirName=''; destination=''
		if len(subdirList) == 0 and os.path.basename(dirName) != 'unclassified' and 'junk' not in str(dirName):
			try:
				base_dirName = os.path.basename(dirName)
				barcode = base_dirName.replace('barcode','')
				sample_name = det_dic[barcode]
				destination = rootDir + '/' + det_dic[barcode]
				new_file = destination+"/"+sample_name+'.fastq' #new_file = destination+"/"+base_dirName+'.fastq'
				os.makedirs(destination, exist_ok=True)
			except:
				continue
			#shutil.move(dirName+"/"+, dest_fpath)
			print([base_dirName, dirName, det_dic[barcode]])
			fileList = [x for x in fileList if x.endswith('.fastq.gz')] #fileList = [x for x in fileList if x.endswith('.fastq')]
			if len(fileList) == 0:
				continue
			n = 0
			bar = progressbar.ProgressBar(maxval=len(fileList), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
			bar.start()
			for filex in fileList:
				os.system('gunzip -c %s >> %s' %(str(dirName)+"/"+filex, new_file)) #base_dirName+'.fastq'))
				n += 1
				bar.update(n)
			bar.finish()
		elif len(subdirList) == 0 and os.path.basename(dirName) == 'x_unclassified':
			fileList = [x for x in fileList if x.endswith('.fastq.gz')]
			destination = rootDir + '/' + 'unclassified_q'
			new_file = destination+'/unclassified.fastq'
			os.makedirs(destination, exist_ok=True)
			for filex in fileList:
				os.system('gunzip -c %s >> %s' %(str(dirName)+"/"+filex, new_file))
	
def return_protocol_dir(protocol):
	schemes = {'artic_v3':'/home/banthony/software/fieldbioinformatics/test-data/primer-schemes',
			  'artic_v4.1':'/home/banthony/software/fieldbioinformatics/artic_v4.1/primer-schemes',
			  'midnight':'/home/banthony/software/fieldbioinformatics/midnight/primer-schemes',
			  'snap':'/home/banthony/software/fieldbioinformatics/snap/primer-schemes',
			  'entebbe':'/home/banthony/software/fieldbioinformatics/entebbe/primer-schemes',
			  'qiaseq':'/home/banthony/software/fieldbioinformatics/qiaseq/primer-schemes',
			  'artic_v4C':'/home/banthony/software/fieldbioinformatics/artic_v4.1/primer-schemes',
			  'artic_v4L':'/home/banthony/software/fieldbioinformatics/artic_v4.1/primer-schemes'}
	#'artic_v4C':'/home/banthony/software/fieldbioinformatics/artic_V4C/primer-schemes'
	return schemes[protocol]

def run_artic(rootDir, details):
	a=b=c=d=e=f=g=x = 0
	destin1 = '/home/banthony/scratch/analysis/covid19/B004_11_2'
	det_dic = read_details(details)

	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0 and 'unclassified' not in str(dirName) and 'barcode' not in str(dirName) and 'NegativeCtrl' not in str(dirName):
			
			fileList = [x for x in fileList if x.endswith('.pchop1.fastq')]
			if len(fileList) == 0:
				continue
			for filex in fileList:
				barcode = filex.split('.')[0].replace('barcode','')
				destination = destin1 + '/' + '/'.join(det_dic[barcode])
				os.makedirs(destination, exist_ok=True)
				protocol, variant, level = det_dic[barcode]
				artic_out= '_'.join([protocol,variant,level,'B'+barcode])
				primer_scheme = return_protocol_dir[protocol]
				read_file= str(dirName)+'/'+filex
				cmd1 = 'cd %s && artic minion --normalise 0 --threads 23 --scheme-directory %s --read-file %s --medaka --medaka-model r941_prom_high_g360 nCoV-2019/V3 %s'
				os.system(cmd1 %(destination, primer_scheme, read_file, artic_out))
				sys.stderr.write(artic_out+'\t'+primer_scheme+'\n')
				
				destination = destination + '/' + 'plot_amplicons'
				os.makedirs(destination, exist_ok=True)
				primer_dir = primer_scheme +  '/nCoV-2019/V3/'
				ref = [x for x in os.listdir(primer_dir) if x.endswith('.reference.fasta')][0]
				bedfile = primer_dir +  'nCoV-2019.bed'
				cmd2 = 'cd %s && samtools ampliconstats --threads 23 --reference %s --max-amplicon-length 3500 %s ../%s.primertrimmed.rg.sorted.bam --output %s.amplicon.stats'
				cmd3 = 'cd %s && plot-ampliconstats -size 1200,900 mydata %s.amplicon.stats'				
				os.system(cmd2 %(destination, primer_dir+ref, bedfile, artic_out, artic_out))
				os.system(cmd3 %(destination, artic_out))
				
				
def collect_metrics(rootDir, details):
	a=b=c=d=e=f=g=x = 0
	destin1 = '/home/banthony/scratch/analysis/covid19/B004_11_2'
	det_dic = read_details(details)

	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0 and 'unclassified' not in str(dirName) and 'barcode' not in str(dirName) and 'NegativeCtrl' not in str(dirName) and 'artic_v3' not in str(dirName):
			
			fileList = [x for x in fileList if x.endswith('.pchop1.fastq')]
			if len(fileList) == 0:
				continue
			for filex in fileList:
				barcode = filex.split('.')[0].replace('barcode','')
				destination = destin1 + '/' + '/'.join(det_dic[barcode]) + '/' + 'metrics'
				os.makedirs(destination, exist_ok=True)
				protocol, variant, level = det_dic[barcode]
				artic_out= '_'.join([protocol,variant,level,'B'+barcode])
				primer_scheme = return_protocol_dir[protocol]
				primer_dir = primer_scheme +  '/nCoV-2019/V3/'
				ref = [x for x in os.listdir(primer_dir) if x.endswith('.reference.fasta')][0]
				ref2 = "/home/banthony/projects/rrg-ioannisr/banthony/analysis/covid19/genome/bwa/NC_045512.2.fasta"
				bedfile = primer_dir +  'nCoV-2019.bed'
				read_file= str(dirName)+'/'+filex
				cmd1 = 'cd %s && nextclade run --reference %s --input-fasta ../%s.consensus.fasta --output-tsv %s_nextclade_results.tsv --input-dataset "/home/banthony/software/sarscov2analysis_SWIFT/data/sars-cov-2"'
				cmd2 = 'cd %s && pangolin %s/%s.consensus.fasta --verbose --outfile %s_pangolin_consensus.csv'
				cmd3 = 'cd %s && minimap2 -ax map-ont -t 23 %s %s | samtools view -bS - | samtools sort -o %s.sorted.bam -'
				cmd4 = 'cd %s &&  samtools flagstat %s.sorted.bam | grep "mapped (" > mapping_rate'
				os.system(cmd1 %(destination, ref2, artic_out, artic_out))
				os.system(cmd2 %(destination, destination.replace('metrics',''),artic_out, artic_out))
				os.system(cmd3 %(destination, primer_dir+ref, read_file, artic_out))
				os.system(cmd4 %(destination, artic_out))
	
def make_adapter_file(rawfile):
	newfile = open(rawfile+'_edited', 'w')
	a=b=c=d=e=f=g=x = 0
	for linex in open(rawfile, 'r'):
		a += 1
		if len(linex.split()) == 2:
			name1 = linex.split()[0]
			name2 = linex.split()[0].strip('_LEFT').strip('_RIGHT')
			seq = linex.strip().split()[-1]
			newfile.write("            Adapter('%s',\n" %name2)
			newfile.write("                    start_sequence=('%s', '%s')),\n" %(name1,seq))
			#newfile.write("                    start_sequence=('%s', '%s'),\n" %(name+'_RIGHT',seq))
		elif len(linex.split()) == 3:
			b += 1
			name1 = linex.split()[0]+'_'+str(b)
			seq1 = linex.strip().split()[-2]
			seq2 = linex.strip().split()[-1]
			newfile.write("            Adapter('%s',\n" %name1)
			newfile.write("                    start_sequence=('%s', '%s'),\n" %(name1+'_RIGHT',seq1))
			newfile.write("                    end_sequence=('%s', '%s')),\n" %(name1+'_LEFT',seq2))			
	newfile.close()
	
def cigar_to_position(samfile):
	#This code takes a samfile and parses it to extract the cigar for each alignment. It tries to return the position on the reference for the mismatches (X), start of deletions (D) and start of insertion (I).
	#The corresponding number of bases for each of the X,D and I is given in first column, followed by IDX then the reference position
	matches = []
	for line_x in open(samfile, 'r'):
		alignment = line_x.strip().split()
		if line_x.startswith('@') or alignment[1] == 4 or alignment[2] == '*':
			continue
		else:
			cigar = alignment[5]
			extra_fields = alignment[11:]
			MD = [i.strip('MD:Z:') for i in extra_fields if 'MD:Z:' in i][0]
			start = 1
			#for num1, i_or_d, num2, m in re.findall('(\d+)([=SIDX])?(\d+)?([A-Za-z])?', cigar):
			#print(re.findall('(\d+)([=SIDXM])', cigar))
			for num1, match in re.findall('(\d+)([=SIDXM])', cigar): #https://stackoverflow.com/questions/17526851/regular-expression-report-position-by-integer-count-from-string
				if match not in 'IDX':
					start += int(num1)
					#if match in 'S':
					#	start += int(num1)
					#continue
				else:
					print(num1, match, start)
					matches += [[num1, match, start]]
					if num1:
						start += int(num1)
					if match in 'I':
						start -= int(num1)
	newfile = open(samfile+'_mutations', 'w')
	for match_x in matches:
		if int(match_x[0]) == 1:
			newfile.write(str(match_x[2])+'\n')
		else:
			if match_x[1] in 'XD':
				for i in range(0,int(match_x[0])):
					newfile.write(str(match_x[2]+i)+'\n')
			#elif match_x[1] == 'X':
			#	for i in range(1,int(match_x[0])):
			#		newfile.write(str(int(match_x[2])+int(match_x[0]))+'\n')
	newfile.close()
	
def manipulate_summary_seq_file(summary_file, detail_file):
	newfile = open(summary_file+'_edited', 'w')
	det_dic = {}; d=0
	print('Creating a dictionary of barcodes and their meaning\n')
	det_dic = read_details(detail_file)
	print('Done creating dictionary. Now counting lines in file ' + summary_file)
	#newfile.write('\t'.join(['protocol','variant','concentration','pass_filter','length','Phred_score','barcode\n']))
	newfile.write('\t'.join(['sample_id','sex','pass_filter','length','Phred_score','barcode', 'sample_id2', 'timepoint\n']))
	length = int(os.popen('wc -l %s' %(summary_file)).read().split()[0])
	print('\nDone counting lines in file ' + summary_file + '. Total number of lines is ...' + str(length+1) + '... \nNow writing final results to file '+ summary_file+'_edited\n')
	bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(summary_file, 'r'):
		line = line.split()
		barcode = line[25]
		pass_filter = 'PASS' if line[10] == 'TRUE' else 'FALSE'
		if barcode.startswith('unclassified'):
			protocol = 'unclassified\tNO'
			newfile.write('\t'.join([ protocol, pass_filter, line[14], str(round(float(line[15]),2)), 'NO', 'NO', 'NO' ])+'\n')
		elif barcode.startswith('barcode') and 'arrangement' not in barcode:
			barcode = barcode.replace('barcode','B')
			try:
				protocol = det_dic[barcode.replace('B','')] #Be very careful, make sure samples in different plates do not share barcodes
			except:
				protocol = 'unclassified\tNO'
				newfile.write('\t'.join([ protocol, pass_filter, line[14], str(round(float(line[15]),2)), barcode, 'NO', 'NO' ])+'\n')
				continue
			sex = protocol[1]; sample = protocol[0]
			embryo = re.findall(r'C\d+(\w){1}',sample)[0]
			timepoint = sample.rstrip(embryo)+'H'
			sample_id2 = timepoint+barcode+sex
			protocol = '\t'.join(protocol)
			newfile.write('\t'.join([ protocol, pass_filter, line[14], str(round(float(line[15]),2)), barcode, sample_id2, timepoint ])+'\n')
		d+=1		
		bar.update(d)
	bar.finish()
	newfile.close()
	
def manipulate_summary_seq_file2(summary_file, detail_file):	
	det_dic = {}; d=0
	print('Creating a dictionary of barcodes and their meaning\n')
	det_dic = read_details2(detail_file)
	if 'B002_05_6_P1' in summary_file:
		det_dic = det_dic['Plate_1']
	elif 'B002_05_6_P2' in summary_file:
		det_dic = det_dic['Plate_2']
	elif 'B002_05_6_P3' in summary_file:
		det_dic = det_dic['Plate_3']
	else:
		print('Cannot properly create the dictionary, please doublecheck..')
		sys.exit()
	
	print('Done creating dictionary. Now counting lines in file ' + summary_file)
	newfile = open(summary_file+'_edited', 'w')
	newfile.write('\t'.join(['sample_id','sex','pass_filter','length','Phred_score','barcode', 'sample_id2', 'timepoint\n']))
	length = int(os.popen('wc -l %s' %(summary_file)).read().split()[0])
	print('\nDone counting lines in file ' + summary_file + '. Total number of lines is ...' + str(length+1) + '... \nNow writing final results to file '+ summary_file+'_edited\n')
	bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(summary_file, 'r'):
		line = line.split()
		barcode = line[25]
		pass_filter = 'PASS' if line[10] == 'TRUE' else 'FAIL'
		if barcode.startswith('unclassified'):
			protocol = 'unclassified\tNO'
			newfile.write('\t'.join([ protocol, pass_filter, line[14], str(round(float(line[15]),2)), 'NO', 'NO', 'NO' ])+'\n')
		elif barcode.startswith('barcode') and 'arrangement' not in barcode:
			barcode = barcode.replace('barcode','B')
			try:
				protocol = det_dic[barcode.replace('B','')] #Be very careful, make sure samples in different plates do not share barcodes
				sex = protocol[1]; sample = protocol[0]
				embryo = re.findall(r'C\d+(\w){1}',sample)[0]
				timepoint = sample.rstrip(embryo)+'H'
				sample_id2 = timepoint+barcode+sex
				protocol = '\t'.join(protocol)
				newfile.write('\t'.join([ protocol, pass_filter, line[14], str(round(float(line[15]),2)), barcode, sample_id2, timepoint ])+'\n')
			except:
				protocol = 'unclassified\tNO'
				newfile.write('\t'.join([ protocol, pass_filter, line[14], str(round(float(line[15]),2)), barcode, 'NO', 'NO' ])+'\n')
				continue
		d+=1		
		bar.update(d)
	bar.finish()
	newfile.close()
	
def entebbe_primers_to_bed(detail_file):
	newfile1 = open(detail_file.replace(os.path.basename(detail_file),'nCoV-2019.bed'), 'w')
	newfile2 = open(detail_file.replace(os.path.basename(detail_file),'nCoV-2019.scheme.bed'), 'w')
	a=b=c=d=e=f=g=h=i=j=k =x= 0
	primer_dic = {}
	for line in open(detail_file, 'r'):
		if line.startswith('Amplicon'):
			continue
		number, primer, sequence, position = line.strip().split()
		if primer.startswith("A") or primer == "nCoV-2019_1_LEFT":
			pool = "nCoV-2019_1"
			if primer.startswith("AF") or 'LEFT' in primer:
				orientation = "_LEFT"
				sign = "+"
			elif primer.startswith("AR"):
				orientation = "_RIGHT"
				sign = "-"
		elif primer.startswith("B") or primer == "nCoV-2019_98_RIGHT":
			pool = "nCoV-2019_2"
			if primer.startswith("BF"):
				orientation = "_LEFT"
				sign = "+"
			elif primer.startswith("BR") or 'RIGHT' in primer:
				orientation = "_RIGHT"
				sign = "-"
		a += 1
		primer_x = "nCoV-2019"+"_"+number+orientation
		to_write = '\t'.join(["NC_045512.2",position,str(int(position)+len(sequence)),"nCoV-2019"+"_"+number+orientation,pool,sign,sequence])+'\n'
		if primer_x not in primer_dic.keys():
			primer_dic.setdefault(primer_x,[to_write])
		else:
			primer_dic[primer_x].append(to_write)
		#newfile1.write('\t'.join(["NC_045512.2",position,str(int(position)+len(sequence)),"nCoV-2019-"+str(a)+"_"+number+orientation,pool,sign,sequence])+'\n')
		#newfile2.write('\t'.join(["NC_045512.2",position,str(int(position)+len(sequence)),"nCoV-2019-"+str(a)+"_"+number+orientation,pool])+'\n')
	for primer,bed in primer_dic.items():
		b = 0
		for line_x in bed:
			ref,start,end,Pid,pool,sign,sequence = line_x.split('\t')
			if b == 0:
				newfile1.write(line_x)
				newfile2.write(line_x.replace('\t'+sign+'\t'+sequence,'\n'))
			else:
				line_x = line_x.replace(Pid,Pid+'_alt'+str(b))
				newfile1.write(line_x)
				newfile2.write(line_x.replace('\t'+sign+'\t'+sequence,'\n'))
			b += 1
	newfile1.close(); newfile2.close()
	
def find_pattern_x(filex):
	match = re.findall('(-\d+_)', cigar)
	
def walk_combine_fastq(rootDir): #I wrote this mainly for snap variant as there was a huge number of unmapped reads for Low and medium
	ref='/home/banthony/software/fieldbioinformatics/entebbe/primer-schemes/nCoV-2019/V3/nCoV-2019.reference.fasta'
	for variant in os.listdir(rootDir):
		print(variant)
		protocol = os.path.basename(rootDir)
		#cmd = 'cd %s && cat High/*pchop1.fastq Low/*pchop1.fastq Medium/*pchop1.fastq > ../%s.%s.all.pchop1.fastq'
		cmd1 = "minimap2 -ax map-ont -t 23 %s %s/Low/*pchop1.fastq | samtools view -S -F 4 - | awk '{print(\">\"$1\"\\n\"$10\"\\n+\\n\"$11)}' > %s.mapped.pchop1.fq" %(ref,rootDir+'/'+variant,rootDir+'/'+variant+'/Low/'+protocol+variant+'Low')
		cmd2 = "minimap2 -ax map-ont -t 23 %s %s/Medium/*pchop1.fastq | samtools view -S -F 4 - | awk '{print(\">\"$1\"\\n\"$10\"\\n+\\n\"$11)}' > %s.mapped.pchop1.fq" %(ref,rootDir+'/'+variant,rootDir+'/'+variant+'/Medium/'+protocol+variant+'Medium')
		cmd3 = 'cd %s && cat High/*pchop1.fastq Low/*mapped.pchop1.fq Medium/*mapped.pchop1.fq > ../%s.%s.all.pchop1.fastq'
		os.system(cmd1)
		os.system(cmd2)
		os.system(cmd3 %(rootDir+'/'+variant, protocol, variant))
		
def return_single_lineseq(file):
	seq = []
	for line in open(file,'r'):
		if not line.startswith('>'):
			seq.append(line.strip())
	return ''.join(seq)
def return_barcode_details(barcode, details):
	det_dic = read_details(details)
	return '.'.join(det_dic[barcode])
def return_average(listx):
	import statistics
	d = 0; total_rds_ampn=[]
	for i in listx:
		if not i.isalpha() and not len(i.split('.')) > 3:
			try:
				total_rds_ampn += [int(i)+0.1]; d += 1
			except:
				total_rds_ampn += [float(i)+0.1]; d += 1
	try:
		stdevn = round(statistics.stdev(total_rds_ampn),3)
		mean_x = round(statistics.mean(total_rds_ampn),3)
		sd_1 = len([x for x in total_rds_ampn if stdevn/x >= 1])
		sd_2 = len([x for x in total_rds_ampn if stdevn/x >= 2])
		sd_3 = len([x for x in total_rds_ampn if stdevn/x >= 3])
		sd_5 = len([x for x in total_rds_ampn if stdevn/x >= 5])
	except:
		mean_x, stdevn, sd_1, sd_2, sd_3, sd_5 = ['.','.','.','.','.','.']
	try:
		CoV = round(stdevn/mean_x,3)
	except:
		CoV = '.'
	return [mean_x, stdevn, CoV, len(total_rds_ampn), sd_1, sd_2, sd_3, sd_5]

def summarise_covseq(listx):
	to_return = []; a=b=c=d=e=f=g=x = 0
	for filex in listx:
		'''	
		artic_out1 = filex.replace('all.pchop1.fastq','')
		protocol = filex.split('.')[0]
		variant = filex.split('.')[-4]
		if 'artic_v4.1' in artic_out1:
			protocol = 'artic_v4.1'
		filex = rootDir+'/'+filex
		'''
		filex = str(dirName)+'/'+filex
		total_number_reads = round(int(os.popen('wc -l %s' %(filex)).read().split()[0])/4,0)
		primer_scheme = return_protocol_dir(protocol)
		primer_dir = primer_scheme +  '/nCoV-2019/V3/'
		ref = [x for x in os.listdir(primer_dir) if x.endswith('.reference.fasta')][0]
		bedfile = primer_dir + 'nCoV-2019.bed'
		ref_length = 29903
		for i in samples:
			if i < total_number_reads:
				b += 1
		print(artic_out1+' to make rounds = ' + str(b+1)+'\n')
		for i in range(b+1):
			#os.system('cd %s && rm -r *' %destination)
			os.makedirs(destination+'plot_amplicons', exist_ok=True)
			if c < b:
				reads_to_sample = samples[i]
			else:
				reads_to_sample = total_number_reads
			c += 1; print('\n'+artic_out1+' round ' + str(c) + ' of '+ str(b+1)+' \t... Sampling ' + str(reads_to_sample) + ' reads.\n')
			cmd1 = 'cd %s && seqtk sample %s %s > R1_sample.fastq'
			cmd2 = 'cd %s && artic minion --normalise 0 --threads 23 --scheme-directory %s --read-file %s --medaka --medaka-model r941_prom_high_g360 nCoV-2019/V3 %s'
			cmd3 = 'cd %s && samtools ampliconstats --threads 23 --reference %s --max-amplicon-length 3500 %s ../%s.primertrimmed.rg.sorted.bam --output %s.amplicon.stats'
			cmd4 = "cd %s && awk 'NR%4==2{a+=length($0)}END{print(a)}' %s > total_bases"
			cmd5 = "cd %s && samtools stats %s.sorted.bam | awk 'FS=\"\t\"{if($2==\"reads mapped:\"){a=$3} else if($2==\"bases mapped (cigar):\"){b=$3}}END{print(a,b)}' > samtools_stats"
			cmd6 = 'cd %s && pangolin %s/%s.consensus.fasta --verbose --outfile %s_pangolin_consensus.csv'
			read_file = 'R1_sample.fastq'
			artic_out = artic_out1.strip(']') + str(reads_to_sample).split('.')[0]
			os.system(cmd1 %(destination,filex,reads_to_sample))
			os.system(cmd2 %(destination, primer_scheme, read_file, artic_out))
			os.system(cmd3 %(destination+'plot_amplicons/', primer_dir+ref, bedfile, artic_out, artic_out))
			#os.system(cmd4 %(destin1,read_file))
			os.system(cmd5 %(destin1,artic_out))
			#os.system(cmd6 %(destination, destination,artic_out, artic_out))
			
			consensus_file = destination+artic_out+".consensus.fasta"
			if os.path.isfile(consensus_file):
				consensus = return_single_lineseq(consensus_file)
				consensus_Ns = consensus.count('N')
				gaps = re.findall(r'(N{5,})', consensus)
				No_gaps = 0; length_gaps = 0
				if len(gaps) != 0:
					for gap in gaps:
						No_gaps += 1; length_gaps += len(gap)
				try:
					ave_len_gap = round(length_gaps/No_gaps,0)
					primersitereport = destin1+'/'+artic_out+'.primersitereport.txt'
					primers_with_mismatch = os.popen('wc -l %s' %(primersitereport)).read().split()[0]
				except:
					consensus_Ns,No_gaps,ave_len_gap,primers_with_mismatch = ['.','.','.','.']
			else:
				consensus_Ns=No_gaps=ave_len_gap=primers_with_mismatch = '.'
			
			ampliconstats_file = amplicon_stats = destination+'/'+'plot_amplicons/'+artic_out+'.amplicon.stats'
			if os.path.isfile(ampliconstats_file):				
				rds_per_ampn = os.popen('grep FREADS %s | tail -1' %(amplicon_stats)).read().split()
				av_rds_per_ampn, stdevn_rds_per_ampn, CoV_rds_per_ampn, No_primers, a1,b1,c1,d1 = return_average(rds_per_ampn)
				bases_per_ampn = os.popen('grep FDEPTH %s | tail -1' %(amplicon_stats)).read().split()
				av_bases_per_ampn, stdevn_bases_per_ampn, CoV_bases_per_ampn, No_primers, a2,b2,c2,d2 = return_average(bases_per_ampn)
				per_per_ampn = os.popen('grep FRPERC %s | tail -1' %(amplicon_stats)).read().split()
				av_per_per_ampn, stdevn_per_per_ampn, CoV_per_per_ampn, No_primers, a3,b3,c3,d3 = return_average(per_per_ampn)
			else:
				av_rds_per_ampn, stdevn_rds_per_ampn, CoV_rds_per_ampn, No_primers, a1,b1,c1,d1 = ['.','.','.','.','.','.','.','.']
				av_bases_per_ampn, stdevn_bases_per_ampn, CoV_bases_per_ampn, No_primers, a2,b2,c2,d2 = ['.','.','.','.','.','.','.','.']
				av_per_per_ampn, stdevn_per_per_ampn, CoV_per_per_ampn, No_primers, a3,b3,c3,d3 = ['.','.','.','.','.','.','.','.']
						
			#total_bases = os.popen('head -1 %s' %(destination+'total_bases')).read().split()[0]
			total_bases = 0
			for i in open(destination+read_file, 'r'):
				e += 1
				if e%4==2:
					total_bases += len(i)
			reads_mapped, bases_mapped = os.popen('head -1 %s' %(destination+'samtools_stats')).read().split()
			
			lineage_file = destination+'/'+artic_out+'_pangolin_consensus.csv'
			if os.path.isfile(lineage_file):
				lineage = os.popen('tail -1 %s' %(lineage_file)).read().split(',')[1]
			else:
				lineage = '.'
			to_write = '\t'.join([protocol,variant,str(reads_to_sample),str(total_bases),reads_mapped,bases_mapped,str(consensus_Ns),str(No_gaps),str(ave_len_gap),str(No_primers),str(primers_with_mismatch),
								  str(av_rds_per_ampn), str(stdevn_rds_per_ampn), str(CoV_rds_per_ampn), str(a1),str(b1),str(c1),str(d1),
								  str(av_bases_per_ampn), str(stdevn_bases_per_ampn), str(CoV_bases_per_ampn), str(a2),str(b2),str(c2),str(d2),
								  str(av_per_per_ampn), str(stdevn_per_per_ampn), str(CoV_per_per_ampn), str(a3),str(b3),str(c3),str(d3),lineage
								 ])
			to_return.append(to_write+'\n')
		b=c=0			
	return to_return

def summarise_snap(listx): #Although I wrote this script, in the end, I am using the artic pipeline but I adjust 2 files
	#File 1: /home/banthony/.local/lib/python3.7/site-packages/artic-1.2.1-py3.7.egg/artic/minion_FOR_SNAP.py line ~146
	#File 2: /home/banthony/.local/lib/python3.7/site-packages/artic-1.2.1-py3.7.egg/artic/make_depth_mask.py lines 58 and 59
	to_return = []; a=b=c=d=e=f=g=x = 0
	for filex in fileList:
		'''
		artic_out1 = filex.replace('all.pchop1.fastq','')
		protocol = filex.split('.')[0]
		variant = filex.split('.')[-4]		
		if 'artic_v4.1' in artic_out1:
			protocol = 'artic_v4.1'
		filex = rootDir+'/'+filex
		'''
		filex = str(dirName)+'/'+filex
		total_number_reads = round(int(os.popen('wc -l %s' %(filex)).read().split()[0])/4,0)
		primer_scheme = return_protocol_dir(protocol)
		primer_dir = primer_scheme +  '/nCoV-2019/V3/'
		ref = [x for x in os.listdir(primer_dir) if x.endswith('.reference.fasta')][0]; ref = primer_dir + ref
		bedfile = primer_dir + 'nCoV-2019.bed'
		for i in samples:
			if i < total_number_reads:
				b += 1
		print(artic_out1+' to make rounds = ' + str(b+1)+'\n')
		for i in range(b+1):
			#os.system('cd %s && rm -r *' %destination)
			os.makedirs(destination+'plot_amplicons', exist_ok=True)
			if c < b:
				reads_to_sample = samples[i]
			else:
				reads_to_sample = total_number_reads
			c += 1; print(artic_out1+' round ' + str(c) + ' of '+ str(b+1)+' \t... Sampling ' + str(reads_to_sample) + ' reads.\n')
			read_file = 'R1_sample.fastq'
			artic_out = artic_out1 + str(reads_to_sample)
			cmd1 = 'cd %s && seqtk sample %s %s > R1_sample.fastq' %(destination,filex,reads_to_sample)
			cmd2 = "cd %s && minimap2 -a -x map-ont -t 23 %s %s | samtools view -bS -F 4 - | samtools sort -o %s.sorted.bam -" %(destination, ref, read_file, artic_out)
			cmd3 = "cd %s && samtools ampliconclip -f %s.clipStats --hard-clip --both-ends --threads 23 -o %s.primclip.sorted.bam --strand -b %s %s" %(destination,artic_out, artic_out, bedfile, artic_out+'.sorted.bam')
			cmd4 = "cd %s && samtools index %s.primclip.sorted.bam" %(destination, artic_out)
			cmd5 = "cd %s && medaka consensus --model r941_prom_high_g360 --threads 2 --chunk_len 300 --chunk_ovlp 100 %s.primclip.sorted.bam %s.hdf" %(destination, artic_out,artic_out)
			cmd6 = "cd %s && medaka variant %s %s.hdf %s.merged.vcf" %(destination, ref, artic_out, artic_out)
			cmd7 = "cd %s && medaka tools annotate --pad 25 %s.merged.vcf %s %s.primclip.sorted.bam tmp.medaka-annotate.vcf" %(destination, artic_out, ref, artic_out)
			cmd8 = "cd %s && mv tmp.medaka-annotate.vcf %s.merged2.vcf" %(destination, artic_out)
			cmd9 = "cd %s && artic_vcf_filter --medaka %s.merged2.vcf %s.pass.vcf %s.fail.vcf" %(destination, artic_out, artic_out, artic_out)
			cmd10 = "cd %s && bgzip -f %s.pass.vcf" %(destination, artic_out)
			cmd11 = "cd %s && tabix -p vcf %s.pass.vcf.gz" %(destination, artic_out)
			cmd12 = "cd %s && samtools view -h --output-fmt SAM %s.primclip.sorted.bam > %s.primclip.sorted.sam" %(destination, artic_out, artic_out)
			cmd13 = "cd %s && samtools view -bS -F 4 %s.primertrimmed.rg.sorted.sam > %s.primertrimmed.rg.sorted.bam" %(destination, artic_out, artic_out)
			cmd14 = "cd %s && samtools index %s.primertrimmed.rg.sorted.bam" %(destination, artic_out)
			cmd15 = "cd %s && artic_make_depth_mask %s %s.primertrimmed.rg.sorted.bam %s.coverage_mask.txt" %(destination, ref, artic_out, artic_out)
			cmd16 = "cd %s && artic_mask %s %s.coverage_mask.txt %s.fail.vcf %s.preconsensus.fasta" %(destination, ref, artic_out, artic_out, artic_out)
			cmd17 = "cd %s && bcftools consensus -f %s.preconsensus.fasta %s.pass.vcf.gz -m %s.coverage_mask.txt -o %s.consensus.fasta" %(destination, artic_out, artic_out, artic_out, artic_out)
			cmd18 = 'cd %s && samtools ampliconstats --threads 23 --reference %s --max-amplicon-length 3500 %s ../%s.primertrimmed.rg.sorted.bam --output %s.amplicon.stats' %(destination+'plot_amplicons/', primer_dir+ref, bedfile, artic_out, artic_out)
			cmd19 = "cd %s && samtools stats %s.sorted.bam | awk 'FS=\"\t\"{if($2==\"reads mapped:\"){a=$3} else if($2==\"bases mapped (cigar):\"){b=$3}}END{print(a,b)}' > samtools_stats" %(destination,artic_out)
			cmd = [cmd1,cmd2,cmd3,cmd4,cmd5,cmd6,cmd7,cmd8,cmd9,cmd10,cmd11,cmd12,cmd13,cmd14,cmd15,cmd16,cmd17,cmd18,cmd19]
			for i in range(19):
				os.system(cmd[i])
				if i == 11:
					samfile = open(destination+'/'+artic_out+'.primertrimmed.rg.sorted.sam','w')
					for line in open(destination+'/'+artic_out+'.primclip.sorted.sam', 'r'):
						f += 1
						if f==2 and line.startswith('@'):
							samfile.write(line+"@RG\tID:1\n")
						elif not line.startswith('@'):
							samfile.write(line.strip()+"\tRG:Z:1\n")
						else:
							samfile.write(line)
					samfile.close(); f = 0
			
			consensus_file = destination+artic_out+".consensus.fasta"
			primersitereport = destin1+'/'+artic_out+'.primersitereport.txt'
			consensus_Ns=No_gaps=ave_len_gap=primers_with_mismatch = '.'
			if os.path.isfile(consensus_file):
				consensus = return_single_lineseq(consensus_file)
				consensus_Ns = consensus.count('N')
				gaps = re.findall(r'(N{5,})', consensus)
				No_gaps = 0; length_gaps = 0
				if len(gaps) != 0:
					for gap in gaps:
						No_gaps += 1; length_gaps += len(gap)
				try:
					ave_len_gap = round(length_gaps/No_gaps,0)					
				except:
					consensus_Ns,No_gaps,ave_len_gap = ['.','.','.']			
			if os.path.isfile(primersitereport):
				try:					
					primers_with_mismatch = os.popen('wc -l %s' %(primersitereport)).read().split()[0]
				except:
					primers_with_mismatch = '.'
			
			ampliconstats_file = amplicon_stats = destination+'/'+'plot_amplicons/'+artic_out+'.amplicon.stats'
			if os.path.isfile(ampliconstats_file):				
				rds_per_ampn = os.popen('grep FREADS %s | tail -1' %(amplicon_stats)).read().split()
				av_rds_per_ampn, stdevn_rds_per_ampn, CoV_rds_per_ampn, No_primers, a1,b1,c1,d1 = return_average(rds_per_ampn)
				bases_per_ampn = os.popen('grep FDEPTH %s | tail -1' %(amplicon_stats)).read().split()
				av_bases_per_ampn, stdevn_bases_per_ampn, CoV_bases_per_ampn, No_primers, a2,b2,c2,d2 = return_average(bases_per_ampn)
				per_per_ampn = os.popen('grep FRPERC %s | tail -1' %(amplicon_stats)).read().split()
				av_per_per_ampn, stdevn_per_per_ampn, CoV_per_per_ampn, No_primers, a3,b3,c3,d3 = return_average(per_per_ampn)
			else:
				av_rds_per_ampn, stdevn_rds_per_ampn, CoV_rds_per_ampn, No_primers, a1,b1,c1,d1 = ['.','.','.','.','.','.','.','.']
				av_bases_per_ampn, stdevn_bases_per_ampn, CoV_bases_per_ampn, No_primers, a2,b2,c2,d2 = ['.','.','.','.','.','.','.','.']
				av_per_per_ampn, stdevn_per_per_ampn, CoV_per_per_ampn, No_primers, a3,b3,c3,d3 = ['.','.','.','.','.','.','.','.']
				
			total_bases = 0
			for i in open(destination+read_file, 'r'):
				e += 1
				if e%4==2:
					total_bases += len(i)
			reads_mapped, bases_mapped = os.popen('head -1 %s' %(destination+'samtools_stats')).read().split()
			
			lineage_file = destination+'/'+artic_out+'_pangolin_consensus.csv'
			if os.path.isfile(lineage_file):
				lineage = os.popen('tail -1 %s' %(lineage_file)).read().split(',')[1]
			else:
				lineage = '.'
			to_write = '\t'.join([protocol,variant,str(reads_to_sample),str(total_bases),reads_mapped,bases_mapped,str(consensus_Ns),str(No_gaps),str(ave_len_gap),str(No_primers),str(primers_with_mismatch),
								  str(av_rds_per_ampn), str(stdevn_rds_per_ampn), str(CoV_rds_per_ampn), str(a1),str(b1),str(c1),str(d1),
								  str(av_bases_per_ampn), str(stdevn_bases_per_ampn), str(CoV_bases_per_ampn), str(a2),str(b2),str(c2),str(d2),
								  str(av_per_per_ampn), str(stdevn_per_per_ampn), str(CoV_per_per_ampn), str(a3),str(b3),str(c3),str(d3),lineage
								 ])
			to_return.append(to_write+'\n')
		b=c=0			
	return to_return

def walk_summarise_covseq(rootDir):
	global dirName,destin1,destination,schemes,artic_out1,protocol,variant,samples,e
	a=b=c=d=e=f=g=x = 0
	destin1=destination = '/home/banthony/scratch/analysis/covid19/B004_twist1/summarise_artic_v3/'
	destin2 = '/home/banthony/scratch/analysis/covid19/B004_twist1/'

	newfile = open(destin2+'summaryresults_artic_v3.txt','w')
	newfile.write('\t'.join(['protocol','variant','No_reads','total_bases','reads_mapped','bases_mapped','consensus_Ns','No_gaps','ave_len_gap','No_primers','No_primers_with_mismatch',
						  'av_rds_per_ampn','stdevn_rds_per_ampn','CoV_rds_per_ampn','rds_1_sd','rds_2_sd','rds_3_sd','rds_5_sd',
						  'av_bases_per_ampn','stdevn_bases_per_ampn','CoV_bases_per_ampn','bs_1_sd','bs_2_sd','bs_3_sd','bs_5_sd',
						  'av_per_per_ampn','stdevn_per_per_ampn','CoV_per_per_ampn','per_1_sd','per_2_sd','per_3_sd','per_5_sd','lineage\n']))
	#fileList = [x for x in os.listdir(rootDir) if x.endswith('all.pchop1.fastq')]
	#samples = [100, 1000, 5000, 10000, 20000, 40000, 80000, 100000, 200000, 400000, 800000, 1000000, 1500000, 2000000]
	#to_write = summarise_covseq(fileList)
	#to_write = summarise_covseq(fileList)
	#newfile.write(to_write+'\n')
	#newfile.close()
	
	samples = [100, 1000, 5000, 10000, 20000, 40000, 80000, 100000, 200000, 400000, 800000, 1000000]
	for dirName, subdirList, fileList in os.walk(rootDir):
		fileList1 = [x for x in os.listdir(dirName) if x.endswith('.subreads.fastq')]
		if len(fileList1) > 0 and len(subdirList) == 0:
		#if len(subdirList) == 0 and not 'artic_v4C' in str(dirName) and not 'artic_v4L' in str(dirName) and not 'NegativeCtrl' in str(dirName) and not 'snapX' in str(dirName):
			#if 'High' in str(dirName) or 'Medium' in str(dirName) or 'Low' in str(dirName):
			if 'Panh' in str(dirName) or 'ONT' in str(dirName) or 'OM' in str(dirName) or 'Wat' in str(dirName) or 'IND' in str(dirName):				
				#fileList1 = [x for x in os.listdir(dirName) if x.endswith('.subreads.fastq')]
				if 'snapX' in str(dirName) and not 'High' in str(dirName):
					fileList1 = [x for x in os.listdir(dirName) if x.endswith('mapped.pchop1.fq')]
				for filey in fileList1:
					#barcode = filex.replace('.pchop1.fastq','').strip('barcode')
					#artic_out1 = return_barcode_details(barcode,'/home/banthony/projects/rrg-ioannisr/banthony/reads/nanopore/covid19/B004/B004_twist1/barcode_assignment')
					artic_out1 =  filey.split('.')[0]+'.' #artic_out1 = '.'.join(str(dirName).strip('/').split('/')[-5:-2])+'.' #artic_out1 = '.'.join(str(dirName).strip('/').split('/')[-3:])+'.' #; print([dirName,artic_out1]) #; print([dirName,artic_out1])
					protocol = 'artic_v4.1' #if 'ONT' in filey else 'Panhandle' #protocol = artic_out1.split('.')[0]
					variant = 'OM:'+artic_out1.strip('.') #'.'.join(artic_out1.split('.')[-3:-1]) #'.'.join(str(artic_out1).split('.')[-3:]) #variant = '.'.join(str(dirName).strip('/').split('/')[-2:]) 
					if 'artic_v4.1' in artic_out1:
						protocol = 'artic_v4.1'				
					fileList2 = [filey]
					to_write = summarise_covseq(fileList2)
					#to_write = summarise_snap(fileList2) #Forget about this line, for SNAP I am using ARTIC pipeline but I change the minion.py and make_depth_mask.py files.
					newfile.write(''.join(to_write))			
	newfile.close()
	
def walk_getLineage_summaries(rootDir):
	a=b=c=d=e=f=g=x = 0
	destination = '/home/banthony/scratch/analysis/covid19/B004_9_4/lineages/'

	newfile = open(destination+'lineage_summaryresult_snap.txt','w')
	newfile.write('\t'.join(['protocol','variant','No_reads','total_skipped','not_correctly_paired','supplementary','skipped_other_rsns','skipped_others',
							 'reads_mapped','mapping_rate','ont_primertrimmed_mapped','otriming_rate','sars_primertrimmed_mapped','ptriming_rate','lineage','conflicts','ambiguity','level\n']))
	
	for dirName, subdirList, fileList in os.walk(rootDir):
		if 'summarise_' in str(dirName) and not 'plot_amplicons' in str(dirName):
			fileList = [x for x in os.listdir(dirName) if x.endswith('.consensus.fasta')]
			for filex in fileList:
				artic_out = filex.replace('.consensus.fasta','')
				protocol = filex.split('.')[0]
				variant = filex.split('.')[1]  #'-'.join(filex.split('.')[1:3])
				level = filex.split('.')[2]
				#try:
				sampled_reads = re.findall('\.(\d{3,})\.', filex)[0] #filex.split('.')[3]
				#if 'artic_v4.1' in filex:
				#	protocol = 'artic_v4.1'
				#	variant = filex.split('.')[2]
				#	sampled_reads = filex.split('.')[3]
				filex = dirName+'/'+filex
				
				consensus_file = filex
				os.system(" cd %s && rm -r pangolin_consensus_%s.csv" %(destination,artic_out))				
				cmd1 = 'cd %s && pangolin %s --verbose --outfile pangolin_consensus_%s.csv'
				os.system(cmd1 %(destination, consensus_file, artic_out))
				lineage_file = destination+'/'+'pangolin_consensus_'+artic_out+'.csv'
				lineage,conflicts,ambiguity = ['.','.','.']
				if os.path.isfile(lineage_file):
					lineage,conflicts,ambiguity = os.popen('tail -1 %s' %(lineage_file)).read().split(',')[1:4]
					if lineage == 'None':
						lineage,conflicts,ambiguity = ['.','.','.']
					if len(ambiguity) < 1:
						ambiguity = '.'
					else:
						ambiguity = str(round(float(ambiguity),3))
						
				align_report_er = dirName+'/'+artic_out+'.alignreport.er'
				total_skipped=supps=no_pair=other=others2 = '.'
				if os.path.isfile(align_report_er):
					total_skipped=supps=no_pair=other=others2 = 0	
					for linex in open(align_report_er, 'r'):
						if not linex.startswith('['):
							total_skipped += 1
						else:
							continue
						if 'skipped as not correctly paired' in linex:
							no_pair += 1
						elif 'skipped as supplementary' in linex:
							supps += 1
						elif 'skipped' in linex:
							other += 1
						else:
							others2 += 1
							
				cmd2 = "samtools flagstats %s | grep 'mapped ('"
				all_mapped='.'; mapping_rate='.'; sars_primertrimmed_mapped='.'; ptriming_rate='.'; ont_primertrimmed_mapped='.'; otriming_rate = '.'
				p1 = dirName+'/'+artic_out+'.sorted.bam'
				p2 = dirName+'/'+artic_out+'.primertrimmed.rg.sorted.bam'
				p3 = dirName+'/'+artic_out+'.trimmed.rg.sorted.bam'
				if os.path.isfile(p1):
					all_mapped = os.popen(cmd2 %p1).read().split()[0].strip('(%') #this is the first mapping done, the number of reads left here
					print([p1,all_mapped,sampled_reads])
					mapping_rate = round((100*int(all_mapped)/int(sampled_reads)),2)
				if os.path.isfile(p2):
					sars_primertrimmed_mapped = os.popen(cmd2 %p2).read().split()[0].strip('(%')
					ptriming_rate = round((100*int(sars_primertrimmed_mapped)/int(sampled_reads)),2)
				if os.path.isfile(p3):
					ont_primertrimmed_mapped = os.popen(cmd2 %p3).read().split()[0].strip('(%')
					otriming_rate = round((100*int(ont_primertrimmed_mapped)/int(sampled_reads)),2)
			
				to_write = '\t'.join([protocol,variant,sampled_reads,str(total_skipped),str(no_pair),str(supps),str(other),str(others2),
								  all_mapped,str(mapping_rate),ont_primertrimmed_mapped,str(otriming_rate),sars_primertrimmed_mapped,str(ptriming_rate),lineage,conflicts,ambiguity,level])
				newfile.write(to_write+'\n')
	newfile.close()
	
def walk_summarise_snap(rootDir):#Although I wrote this script, in the end, I am using the artic pipeline but I adjust 2 files
	#File 1: /home/banthony/.local/lib/python3.7/site-packages/artic-1.2.1-py3.7.egg/artic/minion_FOR_SNAP.py line ~146
	#File 2: /home/banthony/.local/lib/python3.7/site-packages/artic-1.2.1-py3.7.egg/artic/make_depth_mask.py lines 58 and 59
	a=b=c=d=e=f=g=x = 0
	destin1=destination = '/home/banthony/scratch/analysis/covid19/B004_11_2/summarise_snap_1/'
	destin2 = '/home/banthony/scratch/analysis/covid19/B004_11_2/'

	newfile = open(destin2+'summaryresults_snap_1.txt','w')
	newfile.write('\t'.join(['protocol','variant','No_reads','total_bases','reads_mapped','bases_mapped','consensus_Ns','No_gaps','ave_len_gap','No_primers','No_primers_with_mismatch',
						  'av_rds_per_ampn','stdevn_rds_per_ampn','CoV_rds_per_ampn','rds_1_sd','rds_2_sd','rds_3_sd','rds_5_sd',
						  'av_bases_per_ampn','stdevn_bases_per_ampn','CoV_bases_per_ampn','bs_1_sd','bs_2_sd','bs_3_sd','bs_5_sd',
						  'av_per_per_ampn','stdevn_per_per_ampn','CoV_per_per_ampn','per_1_sd','per_2_sd','per_3_sd','per_5_sd','lineage\n']))
	fileList = [x for x in os.listdir(rootDir) if x.endswith('all.pchop1.fastq')]
	samples = [100, 1000, 5000, 10000, 20000, 40000, 80000, 100000, 200000, 400000, 800000, 1000000, 1500000, 2000000]
	for filex in fileList:
		artic_out1 = filex.replace('all.pchop1.fastq','')
		protocol = filex.split('.')[0]
		variant = filex.split('.')[-4]		
		if 'artic_v4.1' in artic_out1:
			protocol = 'artic_v4.1'
		filex = rootDir+'/'+filex
		total_number_reads = int(os.popen('wc -l %s' %(filex)).read().split()[0])/4
		primer_scheme = return_protocol_dir(protocol)
		primer_dir = primer_scheme +  '/nCoV-2019/V3/'
		ref = [x for x in os.listdir(primer_dir) if x.endswith('.reference.fasta')][0]; ref = primer_dir + ref
		bedfile = primer_dir + 'nCoV-2019.bed'
		for i in samples:
			if i < total_number_reads:
				b += 1
		print(artic_out1+' to make rounds = ' + str(b+1)+'\n')
		for i in range(b+1):
			#os.system('cd %s && rm -r *' %destination)
			os.makedirs(destination+'plot_amplicons', exist_ok=True)
			if c < b:
				reads_to_sample = samples[i]
			else:
				reads_to_sample = total_number_reads
			c += 1; print(artic_out1+' round ' + str(c) + ' of '+ str(b+1)+' \t... Sampling ' + str(reads_to_sample) + ' reads.\n')
			read_file = 'R1_sample.fastq'
			artic_out = artic_out1 + str(reads_to_sample)
			cmd1 = 'cd %s && seqtk sample %s %s > R1_sample.fastq' %(destination,filex,reads_to_sample)
			cmd2 = "cd %s && minimap2 -a -x map-ont -t 23 %s %s | samtools view -bS -F 4 - | samtools sort -o %s.sorted.bam -" %(destination, ref, read_file, artic_out)
			cmd3 = "cd %s && samtools ampliconclip -f %s.clipStats --hard-clip --both-ends --threads 23 -o %s.primclip.sorted.bam --strand -b %s %s" %(destination,artic_out, artic_out, bedfile, artic_out+'.sorted.bam')
			cmd4 = "cd %s && samtools index %s.primclip.sorted.bam" %(destination, artic_out)
			cmd5 = "cd %s && medaka consensus --model r941_prom_high_g360 --threads 2 --chunk_len 300 --chunk_ovlp 100 %s.primclip.sorted.bam %s.hdf" %(destination, artic_out,artic_out)
			cmd6 = "cd %s && medaka variant %s %s.hdf %s.merged.vcf" %(destination, ref, artic_out, artic_out)
			cmd7 = "cd %s && medaka tools annotate --pad 25 %s.merged.vcf %s %s.primclip.sorted.bam tmp.medaka-annotate.vcf" %(destination, artic_out, ref, artic_out)
			cmd8 = "cd %s && mv tmp.medaka-annotate.vcf %s.merged2.vcf" %(destination, artic_out)
			cmd9 = "cd %s && artic_vcf_filter --medaka %s.merged2.vcf %s.pass.vcf %s.fail.vcf" %(destination, artic_out, artic_out, artic_out)
			cmd10 = "cd %s && bgzip -f %s.pass.vcf" %(destination, artic_out)
			cmd11 = "cd %s && tabix -p vcf %s.pass.vcf.gz" %(destination, artic_out)
			cmd12 = "cd %s && samtools view -h --output-fmt SAM %s.primclip.sorted.bam > %s.primclip.sorted.sam" %(destination, artic_out, artic_out)
			cmd13 = "cd %s && samtools view -bS -F 4 %s.primertrimmed.rg.sorted.sam > %s.primertrimmed.rg.sorted.bam" %(destination, artic_out, artic_out)
			cmd14 = "cd %s && samtools index %s.primertrimmed.rg.sorted.bam" %(destination, artic_out)
			cmd15 = "cd %s && artic_make_depth_mask %s %s.primertrimmed.rg.sorted.bam %s.coverage_mask.txt" %(destination, ref, artic_out, artic_out)
			cmd16 = "cd %s && artic_mask %s %s.coverage_mask.txt %s.fail.vcf %s.preconsensus.fasta" %(destination, ref, artic_out, artic_out, artic_out)
			cmd17 = "cd %s && bcftools consensus -f %s.preconsensus.fasta %s.pass.vcf.gz -m %s.coverage_mask.txt -o %s.consensus.fasta" %(destination, artic_out, artic_out, artic_out, artic_out)
			cmd18 = 'cd %s && samtools ampliconstats --threads 23 --reference %s --max-amplicon-length 3500 %s ../%s.primertrimmed.rg.sorted.bam --output %s.amplicon.stats' %(destination+'plot_amplicons/', primer_dir+ref, bedfile, artic_out, artic_out)
			cmd19 = "cd %s && samtools stats %s.sorted.bam | awk 'FS=\"\t\"{if($2==\"reads mapped:\"){a=$3} else if($2==\"bases mapped (cigar):\"){b=$3}}END{print(a,b)}' > samtools_stats" %(destination,artic_out)
			cmd = [cmd1,cmd2,cmd3,cmd4,cmd5,cmd6,cmd7,cmd8,cmd9,cmd10,cmd11,cmd12,cmd13,cmd14,cmd15,cmd16,cmd17,cmd18,cmd19]
			for i in range(19):
				os.system(cmd[i])
				if i == 11:
					samfile = open(destination+'/'+artic_out+'.primertrimmed.rg.sorted.sam','w')
					for line in open(destination+'/'+artic_out+'.primclip.sorted.sam', 'r'):
						f += 1
						if f==2 and line.startswith('@'):
							samfile.write(line+"@RG\tID:1\n")
						elif not line.startswith('@'):
							samfile.write(line.strip()+"\tRG:Z:1\n")
						else:
							samfile.write(line)
					samfile.close(); f = 0
			
			consensus_file = destination+artic_out+".consensus.fasta"
			primersitereport = destin1+'/'+artic_out+'.primersitereport.txt'
			consensus_Ns=No_gaps=ave_len_gap=primers_with_mismatch = '.'
			if os.path.isfile(consensus_file):
				consensus = return_single_lineseq(consensus_file)
				consensus_Ns = consensus.count('N')
				gaps = re.findall(r'(N{5,})', consensus)
				No_gaps = 0; length_gaps = 0
				if len(gaps) != 0:
					for gap in gaps:
						No_gaps += 1; length_gaps += len(gap)
				try:
					ave_len_gap = round(length_gaps/No_gaps,0)					
				except:
					consensus_Ns,No_gaps,ave_len_gap = ['.','.','.']			
			if os.path.isfile(primersitereport):
				try:					
					primers_with_mismatch = os.popen('wc -l %s' %(primersitereport)).read().split()[0]
				except:
					primers_with_mismatch = '.'
			
			ampliconstats_file = amplicon_stats = destination+'/'+'plot_amplicons/'+artic_out+'.amplicon.stats'
			if os.path.isfile(ampliconstats_file):				
				rds_per_ampn = os.popen('grep FREADS %s | tail -1' %(amplicon_stats)).read().split()
				av_rds_per_ampn, stdevn_rds_per_ampn, CoV_rds_per_ampn, No_primers, a1,b1,c1,d1 = return_average(rds_per_ampn)
				bases_per_ampn = os.popen('grep FDEPTH %s | tail -1' %(amplicon_stats)).read().split()
				av_bases_per_ampn, stdevn_bases_per_ampn, CoV_bases_per_ampn, No_primers, a2,b2,c2,d2 = return_average(bases_per_ampn)
				per_per_ampn = os.popen('grep FRPERC %s | tail -1' %(amplicon_stats)).read().split()
				av_per_per_ampn, stdevn_per_per_ampn, CoV_per_per_ampn, No_primers, a3,b3,c3,d3 = return_average(per_per_ampn)
			else:
				av_rds_per_ampn, stdevn_rds_per_ampn, CoV_rds_per_ampn, No_primers, a1,b1,c1,d1 = ['.','.','.','.','.','.','.','.']
				av_bases_per_ampn, stdevn_bases_per_ampn, CoV_bases_per_ampn, No_primers, a2,b2,c2,d2 = ['.','.','.','.','.','.','.','.']
				av_per_per_ampn, stdevn_per_per_ampn, CoV_per_per_ampn, No_primers, a3,b3,c3,d3 = ['.','.','.','.','.','.','.','.']
				
			total_bases = 0
			for i in open(destination+read_file, 'r'):
				e += 1
				if e%4==2:
					total_bases += len(i)
			reads_mapped, bases_mapped = os.popen('head -1 %s' %(destination+'samtools_stats')).read().split()
			
			lineage_file = destination+'/'+artic_out+'_pangolin_consensus.csv'
			if os.path.isfile(lineage_file):
				lineage = os.popen('tail -1 %s' %(lineage_file)).read().split(',')[1]
			else:
				lineage = '.'
			to_write = '\t'.join([protocol,variant,str(reads_to_sample),str(total_bases),reads_mapped,bases_mapped,str(consensus_Ns),str(No_gaps),str(ave_len_gap),str(No_primers),str(primers_with_mismatch),
								  str(av_rds_per_ampn), str(stdevn_rds_per_ampn), str(CoV_rds_per_ampn), str(a1),str(b1),str(c1),str(d1),
								  str(av_bases_per_ampn), str(stdevn_bases_per_ampn), str(CoV_bases_per_ampn), str(a2),str(b2),str(c2),str(d2),
								  str(av_per_per_ampn), str(stdevn_per_per_ampn), str(CoV_per_per_ampn), str(a3),str(b3),str(c3),str(d3),lineage
								 ])
			newfile.write(to_write+'\n')
		b=c=0
	newfile.close()
	
def copy_bams_rename(rootDir):
	lin_dic = {'WT':'B.1','UK':'B.1.1.7','SA':'B.1.351','BRZ':'P.1','IND':'B.1.617.2','OM':'BA.1','omi':'BA.1','OMI':'BA.1'}
	files = os.listdir(rootDir)
	reads = ["100","1000","5000","10000","20000","40000","80000","100000"]
	for filex in files:
		#print(filex)
		protocol = filex.replace('summarise_','')
		primer_scheme = return_protocol_dir(protocol)
		primer_dir = primer_scheme +  '/nCoV-2019/V3/'
		ref = [x for x in os.listdir(primer_dir) if x.endswith('.reference.fasta')][0]; ref = primer_dir + ref
		bedfile = primer_dir + 'nCoV-2019.bed'
		destination = '/home/banthony/scratch/analysis/covid19/B004_11_2/plot_amplicons/'+protocol+'/'
		destination2 = '/home/banthony/scratch/analysis/covid19/B004_11_2/plot_amplicons/amp_'+protocol+'/'
		ls_files = os.listdir(rootDir+filex)
		
		for filey in ls_files:
			if filey.endswith('00.sorted.bam'):
				variant = filey.split('.')[-4].split('_')[0]
				sampled_reads = filey.split('.')[-3]
				if sampled_reads in reads:
					lineage = lin_dic[variant]
					newfile = protocol+'-'+lineage+'-'+sampled_reads+'.bam'
					print('Copying ' + filey)
					shutil.copy(rootDir+filex+'/'+filey,destination+newfile)
					shutil.copy(rootDir+filex+'/'+filey+'.bai',destination+newfile+'.bai')
		
		print('\nMaking plots for ' + filex)
		os.system('cd %s && samtools ampliconstats %s %s/*.bam --reference %s --min-depth 20 --max-amplicon-length 3500 --output %s.amplicon.stats --threads 23' %(destination2,bedfile,destination, ref,protocol))
		os.system('cd %s && plot-ampliconstats -size 1200,900 mydata %s.amplicon.stats' %(destination2,protocol))
		
def copy_bams_rename2(rootDir):
	lin_dic = {'WT':'B.1','UK':'B.1.1.7','SA':'B.1.351','BRZ':'P.1','IND':'B.1.617.2','OM':'BA.1','OMI':'BA.1','BR':'P.1','IN':'B.1.617.2'}
	reads = ["100","1000","5000","10000","20000","40000","80000","100000"]
	for dirName, subdirList, fileList in os.walk(rootDir):
		if 'summarise_' in str(dirName) and not 'plot_amplicons' in str(dirName):
			#print(dirName)
			protocol = 'qiaseq'
			destination = '/home/banthony/scratch/analysis/covid19/B004_13_1/plot_amplicons/'+protocol+'/'
			destination2 = '/home/banthony/scratch/analysis/covid19/B004_13_1/plot_amplicons/amp_'+protocol+'/'

			for filey in fileList:
				if filey.endswith('0.sorted.bam') and not filey.endswith('.0.sorted.bam') and '.H.' in filey:
					protocol,variant,level,sampled_reads = filey.split('.')[:4]
					if sampled_reads in reads:
						lineage = lin_dic[variant]
						newfile = protocol+'-'+lineage+'-'+sampled_reads+'.bam'
						#print('Copying ' + filey)
						#shutil.copy(dirName+'/'+filey,destination+newfile)
						#shutil.copy(dirName+'/'+filey+'.bai',destination+newfile+'.bai')
						
	protocol = 'qiaseq'
	destination = '/home/banthony/scratch/analysis/covid19/B004_13_1/plot_amplicons/'+protocol+'/'
	destination2 = '/home/banthony/scratch/analysis/covid19/B004_13_1/plot_amplicons/amp_'+protocol+'/'
	primer_scheme = return_protocol_dir(protocol)
	primer_dir = primer_scheme +  '/nCoV-2019/V3/'
	ref = [x for x in os.listdir(primer_dir) if x.endswith('.reference.fasta')][0]; ref = primer_dir + ref
	bedfile = primer_dir + 'nCoV-2019.bed'
	print('\nMaking plots for ' + protocol)
	#os.system('cd %s && samtools ampliconstats %s %s/*.bam --reference %s --min-depth 20 --max-amplicon-length 3500 --output %s.amplicon.stats --threads 23' %(destination2,bedfile,destination, ref,protocol))
	os.system('cd %s && plot-ampliconstats -size 1200,900 mydata %s.amplicon.stats' %(destination2,protocol))
	
def copy_bams_rename3(rootDir,rootDir2):
	lin_dic = {'WT':'B.1','UK':'B.1.1.7','SA':'B.1.351','BRZ':'P.1','IND':'B.1.617.2','OM':'BA.1','omi':'BA.1','OMI':'BA.1'}
	files = os.listdir(rootDir)
	reads = ["100","1000","5000","10000","20000","40000","80000","100000"]
	'''
	for filex in files:		
		protocol = filex.split('.')[0] if not 'v4.1' in filex else 'artic_v4.1'
		destination = '/home/banthony/scratch/analysis/covid19/B004_11_2/plot_amplicons/'+protocol+'/'
		destination2 = '/home/banthony/scratch/analysis/covid19/B004_11_2/plot_amplicons/amp_'+protocol+'/'
		ls_files = [filex]
		
		for filey in ls_files:
			if filey.endswith('00.sorted.bam') and 'Wat' not in filey and 'snap' not in filey and 'vLow' not in filey:
				#print(filey)
				variant = filey.split('.')[-5]
				sampled_reads = filey.split('.')[-3]
				if sampled_reads in reads:
					lineage = lin_dic[variant]
					newfile = protocol+'-'+lineage+'-'+sampled_reads+'.bam'
					print('Copying ' + filey)
					shutil.copy(rootDir+'/'+filey,destination+newfile)
					shutil.copy(rootDir+'/'+filey+'.bai',destination+newfile+'.bai')
	
	files = os.listdir(rootDir2)
	for filey in files:
		protocol = 'snap'
		destination = '/home/banthony/scratch/analysis/covid19/B004_11_2/plot_amplicons/'+protocol+'/'
		if filey.endswith('00.sorted.bam') and 'Wat' not in filey and 'vLow' not in filey:
			variant = filey.split('.')[-5]
			sampled_reads = filey.split('.')[-3]
			if sampled_reads in reads:
				lineage = lin_dic[variant]
				newfile = protocol+'-'+lineage+'-'+sampled_reads+'.bam'
				print('Copying ' + filey)
				shutil.copy(rootDir2+'/'+filey,destination+newfile)
				shutil.copy(rootDir2+'/'+filey+'.bai',destination+newfile+'.bai')
	'''			
	for protocol in os.listdir('/home/banthony/scratch/analysis/covid19/B004_11_2/plot_amplicons/'):
		if 'amp_' in protocol or 'analysis' in protocol:
			continue
		primer_scheme = return_protocol_dir(protocol)
		primer_dir = primer_scheme +  '/nCoV-2019/V3/'
		ref = [x for x in os.listdir(primer_dir) if x.endswith('.reference.fasta')][0]; ref = primer_dir + ref
		bedfile = primer_dir + 'nCoV-2019.bed'
		destination = '/home/banthony/scratch/analysis/covid19/B004_11_2/plot_amplicons/'+protocol+'/'
		destination2 = '/home/banthony/scratch/analysis/covid19/B004_11_2/plot_amplicons/amp_'+protocol+'/'
		print('\nMaking plots for ' + protocol)
		os.system('cd %s && samtools ampliconstats %s %s/*.bam --reference %s --min-depth 20 --max-amplicon-length 3500 --output %s.amplicon.stats --threads 23' %(destination2,bedfile,destination, ref,protocol))
		os.system('cd %s && plot-ampliconstats -size 1200,900 mydata %s.amplicon.stats' %(destination2,protocol))

def return_total_reads(rootDir, details):
	det_dic = read_details(details)
	new_dic ={}
	new_dic2 ={}
	for dirName, subdirList, fileList in os.walk(rootDir):
		#if 'midnight' not in str(dirName) or 'entebbe' not in str(dirName) or 'snap' not in str(dirName) or 'artic_v3' not in str(dirName) or 'artic_v4.1' not in str(dirName):
		#	continue
		if len(subdirList) == 0 and os.path.basename(dirName) != 'unclassified_q' and len([x for x in fileList if x.endswith('pchop1.fastq')]) >0:
			print([dirName])
			if not 'snap' in str(dirName):
				#print('Starting..')
				fileList = [x for x in fileList if x.endswith('pchop1.fastq')]
				for filex in fileList:
					length = int(os.popen('wc -l %s' %(str(dirName)+'/'+filex)).read().split()[0])
					reads = round((length)/4,0)
					new_dic.setdefault(filex.split('.')[0].replace('barcode',''), det_dic[filex.split('.')[0].replace('barcode','')]+ [str(reads)+'\n'])
			elif 'snap' in str(dirName):
				if 'High' in str(dirName):
					fileList = [x for x in fileList if x.endswith('pchop1.fastq')]
					for filex in fileList:
						length = int(os.popen('wc -l %s' %(str(dirName)+'/'+filex)).read().split()[0])
						reads = round((length)/4,0)
						new_dic.setdefault(filex.split('.')[0].replace('barcode',''), det_dic[filex.split('.')[0].replace('barcode','')]+ [str(reads)+'\n'])
				elif 'Low' in str(dirName) or 'Medium' in str(dirName):
					fileList1 = [x for x in fileList if x.endswith('pchop1.fastq')]
					fileList2 = [x for x in fileList if x.endswith('mapped.pchop1.fq')][0]
					for filex in fileList1:
						length = int(os.popen('wc -l %s' %(str(dirName)+'/'+fileList2)).read().split()[0])
						reads = round((length)/4,0)
						new_dic.setdefault(filex.split('.')[0].replace('barcode',''), det_dic[filex.split('.')[0].replace('barcode','')]+ [str(reads)+'\n'])
			#print(str(dirName)+' done, moving on...')
	newfile = open(rootDir+'/number_of_pchop_reads.tsv', 'w')
	for barcode,details in new_dic.items():
		newfile.write('\t'.join(details))

def return_read_length_list(filex):
	read_length_list = []
	for tit,det in fastq_to_dict(filex).items():
		read_length = str(len(det.split('\n')[0]))
		read_length_list += [read_length]
	return read_length_list

def return_read_lengths_details(rootDir, details):
	det_dic = read_details(details)
	new_dic = {}
	for dirName, subdirList, fileList in os.walk(rootDir):
		a=1; dirs = ['BR','IN','OM','SA','UK','WT']
		if len(subdirList) == 0 and os.path.basename(dirName) != 'unclassified_q' and len([x for x in fileList if x.endswith('pchop1.fastq')]) > 0:
			#if 'BR' not in str(dirName) or 'IN' not in str(dirName) or 'OM' not in str(dirName) or 'SA' not in str(dirName) or 'UK' not in str(dirName) or 'WT' not in str(dirName):
			if 'LSPQ' in str(dirName) or 'vLow'  in str(dirName) or 'Wat'  in str(dirName):
				continue
			if 'BR' in str(dirName) or 'IN' in str(dirName) or 'OM' in str(dirName) or 'SA' in str(dirName) or 'UK' in str(dirName) or 'WT' in str(dirName):
			#if 'artic_v3' in str(dirName) or 'artic_v4.1' in str(dirName) or 'entebbe' in str(dirName) or 'midnight' in str(dirName) or 'snap' in str(dirName):
				print([dirName])
				fileList = [x for x in fileList if x.endswith('pchop1.fastq')]
				for filex in fileList:
					a+=1
					barcode_details = '\t'.join(det_dic[filex.split('.')[0].replace('barcode','')])
					if barcode_details not in new_dic.keys():
						new_dic.setdefault(barcode_details,return_read_length_list(dirName+'/'+filex))
					elif barcode_details in new_dic.keys():
						new_dic.setdefault(barcode_details+str(a),return_read_length_list(dirName+'/'+filex))
			elif 'snap' in str(dirName):
				if 'High' in str(dirName):
					fileList = [x for x in fileList if x.endswith('pchop1.fastq')]
					for filex in fileList:
						barcode_details = '\t'.join(det_dic[filex.split('.')[0].replace('barcode','')])
						if barcode_details not in new_dic.keys():
							new_dic.setdefault(barcode_details,return_read_length_list(dirName+'/'+filex))
				elif 'Low' in str(dirName) or 'Medium' in str(dirName):
					fileList1 = [x for x in fileList if x.endswith('pchop1.fastq')]
					fileList2 = [x for x in fileList if x.endswith('mapped.pchop1.fq')][0]
					for filex in fileList1:
						barcode_details = '\t'.join(det_dic[filex.split('.')[0].replace('barcode','')])
						if barcode_details not in new_dic.keys():
							new_dic.setdefault(barcode_details,return_read_length_list(dirName+'/'+fileList2))
	newfile = open(rootDir+'/B004_11_3_pass_read_lengths_details.tsv', 'w')
	for barcode,details in new_dic.items():
		for i in details:
			newfile.write(barcode + '\t' + i + '\n')
	newfile.close()
	
def return_read_lengths_details2(rootDir, details):
	det_dic = {}
	for i in open(details, 'r'):
		det_dic.setdefault(i.split()[1],i.split()[0]+'\t'+i.split()[2])
	new_dic = {}
	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0 and os.path.basename(dirName) != 'unclassified_q' and len([x for x in fileList if x.endswith('.all.porechop.pychop.cutadapt.fq')]) > 0:
			print(dirName)
			fileList = [x for x in fileList if x.endswith('.all.porechop.pychop.cutadapt.fq')]
			for filex in fileList:
				sample_id = filex.split('.')[0]
				embryo = re.findall(r'C\d+(\w){1}',sample_id)[0]
				timepoint = sample_id.rstrip(embryo)+'H'
				barcode = det_dic[sample_id].split('\t')[0]; sex = det_dic[sample_id].split('\t')[1]
				sample_id2 = timepoint+barcode+sex
				barcode_details = '\t'.join([sample_id2,sample_id,barcode,sex])
				new_dic.setdefault(barcode_details,return_read_length_list(dirName+'/'+filex))
	print('Done creating dictionary. Now writing final results to file...')			
	newfile = open(rootDir+'/processed_read_lengths', 'w')
	for barcodex,details in new_dic.items():
		for i in details:
			newfile.write(barcodex + '\t' + i + '\n')
	newfile.close()
	
def correct_qiaseq(qiaseq_bed):
	dirc = os.path.dirname(qiaseq_bed)
	nfile1 = open(dirc+'/'+'SARS-CoV-2.scheme.bed', 'w')
	for line in open(qiaseq_bed, 'r'):
		match = re.findall('(-\d+_)', line)
		if match:
			line = line.replace(match[0],'_')
			line = line.replace(line.split()[3],line.split()[3]+'_alt'+match[0].strip('-_'))
			nfile1.write(line)
		else:
			nfile1.write(line)
			
def quick_merge(rootDir):
	dirs = ['artic_v3','artic_v4.1','artic_v4C','artic_v4L','entebbe','midnight','snap']
	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0:
			for i in dirs:
				if i in dirName:
					fileList = [x for x in fileList if x.endswith('.fastq')]
					cmd = 'cat %s >> %s'
					for filex in fileList:
						if os.path.isfile(dirName.replace('fastq_pass2','fastq_pass')+'/'+filex):
							print('Merging...: ' + dirName+'/'+filex + ' and ' + dirName.replace('fastq_pass2','fastq_pass')+'/'+filex)
							os.system(cmd %(dirName+'/'+filex, dirName.replace('fastq_pass2','fastq_pass')+'/'+filex))
						else:
							print('File does not exist at destination: ' + dirName+'/'+filex)
			
def R2C2(rootDir):
	dirs = ['artic_v4C']
	splint_file='/home/banthony/projects/rrg-ioannisr/banthony/analysis/covid19/B004_09_3/lambda_phage_splint.fasta'
	adapters='/home/banthony/projects/rrg-ioannisr/banthony/analysis/covid19/B004_09_3/adapter.fasta'
	config_file='/home/banthony/scratch/analysis/covid19/B004_9_4/config'
	NUM_mat='/home/banthony/software/C3POa_gonk/NUC.4.4.mat'
	for dirName, subdirList, fileList in os.walk(rootDir):
		fileList = [x for x in fileList if x.endswith('.pchop1.fastq')]
		if len(fileList) > 0:			
		#if len(subdirList) == 0:
			#fileList = [x for x in fileList if x.endswith('.pchop1.fastq')]
			for filex in fileList:
				reads_file=dirName+'/'+filex
				print('R2C2 for...: ' + reads_file)
				barcode = filex.strip('.pchop1.fastq')
				pre_process_dir = dirName+'/R2C2_preprocess_'+barcode; os.makedirs(pre_process_dir, exist_ok=True)
				pre_process_output = pre_process_dir+'/lambda_phage_splint/R2C2_raw_reads.fastq'
				process_output = dirName+'/R2C2_process_'+barcode; os.makedirs(process_output, exist_ok=True)
				
				os.system('cd %s && rm -r *' %pre_process_dir)
				cmd1 = 'python /home/banthony/software/C3POa_gonk/C3POa_preprocessing.py -i %s -o %s -q %s -l %s -s %s -c %s' %(reads_file, pre_process_dir, '3', '300', splint_file, config_file)
				#cmd2 = 'python /home/banthony/software/C3POa_gonk/C3POa.py -p %s -m %s -l %s -c %s -r %s -n 1' %(process_output, NUM_mat, '1000', config_file, pre_process_output)
				os.system(cmd1)
				#os.system(cmd2)
				
				
def split_R2C2(fastq):
	length = int(os.popen('wc -l %s' %(fastq)).read().split()[0])
	file1=open(fastq,'r')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	x=0
	newfile = open(fastq.replace('.fastq','.subreads.fastq'),'w')
	while x < length:
		name=file1.readline()
		seq=file1.readline()
		plus=file1.readline()
		qual=file1.readline()
		
		a=0
		splint_pos = name.strip().split('_')[1:]
		pos = len(splint_pos)
		for i in range(len(splint_pos)+1):
			if i == 0:
				left = 0; right = splint_pos[i]
			elif i > 0 and i < len(splint_pos):
				left = splint_pos[i-1]; right = splint_pos[i]
			else:
				left = splint_pos[i-1]; right = len(seq)
			a += 1	
			new_name = name.strip().split('_')[0]+'_'+str(a)+'\n'
			new_seq = seq.strip()[int(left):int(right)]
			new_qual = qual.strip()[int(left):int(right)]
			newfile.write(new_name+new_seq+'\n'+plus+new_qual+'\n')
		x+=4
		bar.update(x)
	bar.finish()
	newfile.close()
	
def walk_split_R2C2(rootDir):
	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0:
			fileList = [x for x in fileList if x.endswith('R2C2_raw_reads.fastq')]
			for filex in fileList:				
				pre_process_output = dirName+'/'+filex
				print(pre_process_output)
				split_R2C2(pre_process_output)
				
def artic_to_custom(artic_bed):
	newfile1 = open(artic_bed+'.custom.tsv', 'w')
	for line in open(artic_bed, 'r'):
		primer_name = line.split()[3]
		sequence = line.strip().split()[6]
		if 'LEFT' in line:
			newfile1.write(primer_name.replace('LEFT','LEFTcus') + '\t' + 'AATGATACGGCGACCACCGAGATCTACAC' + sequence + '\n')
		elif 'RIGHT' in line:
			newfile1.write(primer_name.replace('RIGHT','RIGHTcus') + '\t' + 'AAGCAGTGGTATCAACGCAGAGT' + sequence + '\n')
	newfile1.close()
	
def lampore_assign(read_annotation_dir, rootdir):
	a=0; c=0
	read_annot = {}
	to_write = {}
	annot_files = [x for x in os.listdir(read_annotation_dir) if x.endswith('read_annotation.tsv')]
	print('Collecting barcode-FIP data from read annotation files')
	bar = progressbar.ProgressBar(maxval=len(annot_files), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for filex in annot_files:
		for line in open(read_annotation_dir+'/'+filex, 'r'):
			read_id = line.split()[0]
			barcode = line.split()[-1]
			if not read_id in read_annot.keys():
				read_annot.setdefault(read_id,barcode)
			#else:
			#	print(read_id)
		a += 1
		bar.update(a)
	bar.finish()
	print('Assigning reads to barcode-FIPs')
	for dirName, subdirList, fileList in os.walk(rootdir):
		b=0
		if len(subdirList) == 0:
			fileList = [x for x in fileList if x.endswith('.fastq')]
			if len(fileList) != 0:
				print(dirName) 
			else:
				continue
			bar = progressbar.ProgressBar(maxval=len(fileList), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
			bar.start()
			for filex in fileList:
				for tit,det in fastq_to_dict(dirName+'/'+filex).items():
					if tit.strip() in read_annot.keys():
						barcode_x = read_annot[tit.strip()]
						if barcode_x not in to_write.keys():
							to_write.setdefault(barcode_x, [det])
						else:
							to_write[barcode_x].append(det)
				b += 1
				bar.update(b)
			bar.finish()
	print('Writing final fastq files')
	destin = rootdir+'/barcode_FIP'
	os.makedirs(destin, exist_ok=True)
	newfile1 = open(destin+'/reads_manifest.tsv','w')
	bar = progressbar.ProgressBar(maxval=len(to_write.keys()), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for barcodes,reads in to_write.items():
		newfile = open(destin+'/'+barcodes+'.fastq','w')
		newfile1.write('\t'.join([barcodes.split('_')[0],barcodes.split('_')[1],'x',destin+'/'+barcodes+'.fastq\n']))
		newfile.write(''.join(reads))
		newfile.close()
		c += 1
		bar.update(c)
	bar.finish()
	newfile1.close()
	
def lampore_analysis(rootdir):
	'''
	#Traverse rootdir to get to fastq_pass then assign reads to barcode-FIP
	for dirName, subdirList, fileList in os.walk(rootdir):
		if 'fastq_pass' in subdirList and 'fastq_fail' in subdirList and 'lampore' in subdirList:
			if os.path.isdir(dirName + '/fastq_pass/barcode_FIP'):
				continue
			sample_id = re.findall(r'EPO_(PL0\d+\w?)_', dirName)[0]
			fastq_pass_dir = dirName + '/fastq_pass/'
			read_annotation_dir = dirName + '/lampore/' + os.listdir(dirName + '/lampore')[0] + '/vsearch/'
			lampore_assign(read_annotation_dir, fastq_pass_dir)
			
			#Call flair conda deactivate && module load python/3.10.2 && module load samtools/1.10 && module load bedtools/2.30.0 && export PATH=/home/banthony/software/minimap2-2.24_x64-linux:$PATH
			isoforms='/home/banthony/projects/rrg-ioannisr/banthony/analysis/covid19/genome/chromFa/lampore_ref_ed.fasta'
			reads_manifest=fastq_pass_dir+'/barcode_FIP/reads_manifest.tsv'
			os.system('cd %s && python /home/banthony/software/flair/flair.py quantify --threads 23 --temp_dir tmp -r %s -i %s --output %s' %('/home/banthony/scratch/analysis/covid19/EPO/flair_output',reads_manifest,isoforms, sample_id+'_counts_matrix.tsv'))
			os.system('rm -r /home/banthony/scratch/analysis/covid19/EPO/flair_output/tmp')
	'''
	#evaluate flair and combine all flair results	
	newfile = open('/home/banthony/scratch/analysis/covid19/EPO/pass_flair_output.tsv','w')
	newfile1 = open('/home/banthony/scratch/analysis/covid19/EPO/pass_flair_output_invalid.tsv','w')
	newfile.write('\t'.join(['plate','ids','ACTB','E','N','ORF1a','ORF1ab','result\n']))
	newfile1.write('\t'.join(['plate','ids','ACTB','E','N','ORF1a','ORF1ab','result\n']))
	for file_x in [y for y in os.listdir('/home/banthony/scratch/analysis/covid19/EPO/flair_output') if y.endswith('counts_matrix.tsv')]:
		print(file_x)
		sample_idx = file_x.split('_')[0]
		flair_dic = {}
		for line in open('/home/banthony/scratch/analysis/covid19/EPO/flair_output/'+file_x, 'r'):
			flair_dic.setdefault(line.split()[0],line.split()[1:])
		for i in range(len(flair_dic['ids'])):
			sample_id_x = flair_dic['ids'][i]
			if sample_id_x in ['barcode12_FIP05_x','barcode12_FIP06_x','barcode12_FIP07_x','barcode12_FIP08_x']:
				continue
			ACTB_genomic,ACTB_transcript,E,N,ORF1a,ORF1ab = ['0','0','0','0','0','0']
			try:
				ACTB_genomic = flair_dic['ACTB_genomic'][i].split('.')[0]
				ACTB_transcript = flair_dic['ACTB_transcript'][i].split('.')[0]
				E = flair_dic['E'][i].split('.')[0]
				N = flair_dic['N'][i].split('.')[0]
				ORF1a = flair_dic['ORF1a'][i].split('.')[0]
				ORF1ab = flair_dic['ORF1ab'][i].split('.')[0]
			except:
				k = 1
			if int(ACTB_transcript) > 49 and (int(E) + int(N) + int(ORF1a) + int(ORF1ab)) == 0: #Valid run and test is Negative
				continue
			elif int(ACTB_transcript) < 50 and (int(E) + int(N) + int(ORF1a) + int(ORF1ab)) == 0:
				lamp = 'Invalid'; newfile1.write('\t'.join([sample_idx,sample_id_x,ACTB_transcript,E,N,ORF1a,ORF1ab,lamp+'\n']))
			else:
				#elif int(ACTB_transcript) > 49:
				if int(E) + int(N) + int(ORF1a) >= 50 or int(E) + int(N) + int(ORF1ab) >= 50:
					lamp = 'Positive'; newfile.write('\t'.join([sample_idx,sample_id_x,ACTB_transcript,E,N,ORF1a,ORF1ab,lamp+'\n']))
				elif int(E) + int(N) + int(ORF1a) >= 20 or int(E) + int(N) + int(ORF1ab) >= 20:
					lamp = 'inconclusive'; newfile.write('\t'.join([sample_idx,sample_id_x,ACTB_transcript,E,N,ORF1a,ORF1ab,lamp+'\n']))
				elif int(E) + int(N) + int(ORF1a) < 20 or int(E) + int(N) + int(ORF1ab) < 20:
					lamp = 'Negative'; newfile.write('\t'.join([sample_idx,sample_id_x,ACTB_transcript,E,N,ORF1a,ORF1ab,lamp+'\n']))
	newfile.close(); newfile1.close()
	
def combine_liqa_results(rootdir): #, gene_transcript):
	liqa_comb = {}
	sample = {}
	newfile = open(rootdir+'/MoY_liqa','w')
	newfile.write('\t'.join(['sample_id','GeneName','IsoformName','ReadPerGene_corrected','RelativeAbundance','infor_ratio\n']))
	for dirName, subdirList, fileList in os.walk(rootdir):
		if len(subdirList) == 0 and 'barcode' not in dirName:
			b_dirName = os.path.basename(dirName.rstrip('/'))
			fileList2 = [x for x in fileList if x.endswith('.all.fastq')]
			fileList = [x for x in fileList if x.endswith('porechop.liqa')]
			if len(fileList) > 0:
				for filex in fileList:
					for line in open(dirName+'/'+filex, 'r'):
						if 'PB.CAJ.269.1' in line or 'PB.CAJ.269.2' in line:
							newfile.write(filex.split('.')[0]+'\t'+line)
				#for line in open(gene_transcript, 'r'):
				#	liqa = l
	newfile.close()
				
def edit_gtf(gtf):
	newfile = open(gtf+'_ed','w')
	for line in open(gtf, 'r'):
		if "gene_name" not in line and  'gene_id' in line:
			newfile.write(line.strip() + ' gene_name "%s";' %(re.findall(r'gene_id "(\S+)";',line)[0]) + '\n')
		elif "gene_name" in line and  'gene_id' in line:
			gene_id = re.findall(r'gene_id "(\S+)";',line)[0]
			gene_name = re.findall(r'gene_name "(\S+)";',line)[0]
			if gene_name.startswith('LOC') and not gene_id.startswith('LOC'):
				newfile.write(line.replace(gene_name, gene_id) )
			else:
				newfile.write(line)
		else:
			newfile.write(line)
	newfile.close()
	
def corinne_file_man(corinne_file):
	newfile = open(corinne_file+'_ed','w')
	for line in open(corinne_file, 'r'):
		a = 0
		Lampore_barcode,Genome_completeness,Epo_ID_code_for_Pos_or_inconclusive_samples,Resultat_du_test_rapide,Reultat_PCR_gargarisme,Resultat_PCR_nasal,Postive_RT_PCR = line.split('\t')[:7]
		lamp_result = line.split('\t')[-3]
		MGC_qPCR_results = line.split('\t')[-1]
		#print([Lampore_barcode,Genome_completeness,Epo_ID_code_for_Pos_or_inconclusive_samples,Resultat_du_test_rapide,Reultat_PCR_gargarisme,Resultat_PCR_nasal,Postive_RT_PCR])
		if line.startswith('Lampore_barcode'):
			newfile.write(line.strip()+'\tTest\n')
			continue
		if Genome_completeness == 'NA' :
			Genome_completeness = 0; a += 1
		if float(Genome_completeness) > 0:
			Genome_completeness = 1
		if Resultat_du_test_rapide == 'NA':
			Resultat_du_test_rapide = 0; a += 1
		if Reultat_PCR_gargarisme == 'NA':
			Reultat_PCR_gargarisme = 0; a += 1
		if Resultat_PCR_nasal == 'NA' or Resultat_PCR_nasal == '':
			Resultat_PCR_nasal = 0; a += 1
		if Postive_RT_PCR == 'NA' or Postive_RT_PCR == '':
			Postive_RT_PCR = 0; a += 1
		#print([int(Genome_completeness) , Resultat_du_test_rapide , Reultat_PCR_gargarisme , Resultat_PCR_nasal , Postive_RT_PCR])
		#True positive
		#if MGC_qPCR_results == 'Positive' or int(Genome_completeness) > 0 or (int(Resultat_du_test_rapide) + int(Reultat_PCR_gargarisme) + int(Resultat_PCR_nasal) + int(Postive_RT_PCR) >= 2):
		if MGC_qPCR_results == 'Positive' or int(Genome_completeness) > 0 or int(Resultat_du_test_rapide) > 0 or int(Reultat_PCR_gargarisme) > 0 or int(Resultat_PCR_nasal) > 0 or int(Postive_RT_PCR) > 0:
			if lamp_result == 'Positive':
				newfile.write(line.strip()+'\tTrue positive\n')
			elif lamp_result == 'Negative':
				newfile.write(line.strip()+'\tFalse negative\n')

		elif MGC_qPCR_results == 'Negative' or ((int(Genome_completeness) + int(Resultat_du_test_rapide) + int(Reultat_PCR_gargarisme) + int(Resultat_PCR_nasal) + int(Postive_RT_PCR) == 0) and a < 4):
			if lamp_result == 'Positive':
				newfile.write(line.strip()+'\tFalse positive\n')
			elif lamp_result == 'Negative':
				newfile.write(line.strip()+'\tTrue negative\n')
		#if int(Genome_completeness) + int(Resultat_du_test_rapide) + int(Reultat_PCR_gargarisme) + int(Resultat_PCR_nasal) + int(Postive_RT_PCR) >= 2:
			#newfile.write(line)
	newfile.close()
	
def corinne_file_man2(corinne_file): #Please ignore this code, it needs so much improvement on the True negatives
	newfile = open(corinne_file+'_ed2','w')
	for line in open(corinne_file, 'r'):
		a = 0
		Lampore_barcode,Genome_completeness,Epo_ID_code_for_Pos_or_inconclusive_samples,Resultat_du_test_rapide,Reultat_PCR_gargarisme,Resultat_PCR_nasal,Postive_RT_PCR = line.split('\t')[:7]
		Lampore_barcode1,Genome_completeness1,Epo_ID_code_for_Pos_or_inconclusive_samples1,Resultat_du_test_rapide1,Reultat_PCR_gargarisme1,Resultat_PCR_nasal1,Postive_RT_PCR1 = line.split('\t')[:7]
		lamp_result = line.split('\t')[-3]
		MGC_qPCR_results = line.split('\t')[-1]
		MGC_qPCR_results1 = 1 if MGC_qPCR_results == 'Positive' else 0
		#print([Lampore_barcode,Genome_completeness,Epo_ID_code_for_Pos_or_inconclusive_samples,Resultat_du_test_rapide,Reultat_PCR_gargarisme,Resultat_PCR_nasal,Postive_RT_PCR])
		if line.startswith('Lampore_barcode'):
			newfile.write(line.strip()+'\tTest\n')
			continue		
		if Genome_completeness == 'NA' :
			Genome_completeness = 0; a += 1
		if float(Genome_completeness) > 0:
			Genome_completeness = 1
		if Resultat_du_test_rapide == 'NA':
			Resultat_du_test_rapide = 0; a += 1
		if Reultat_PCR_gargarisme == 'NA':
			Reultat_PCR_gargarisme = 0; a += 1
		if Resultat_PCR_nasal == 'NA' or Resultat_PCR_nasal == '':
			Resultat_PCR_nasal = 0; a += 1
		if Postive_RT_PCR == 'NA' or Postive_RT_PCR == '':
			Postive_RT_PCR = 0; a += 1
		MGC_qPCR_results2 = 0 if MGC_qPCR_results == 'NA' else 1; a += 1
		#True positive
		if MGC_qPCR_results == 'Positive' or int(Genome_completeness) > 0 or int(Resultat_du_test_rapide) > 0 or int(Reultat_PCR_gargarisme) > 0 or int(Resultat_PCR_nasal) > 0 or int(Postive_RT_PCR) > 0:
			if lamp_result == 'Positive':
				newfile.write(line.strip()+'\tTrue positive\n')
			elif lamp_result == 'Negative':
				newfile.write(line.strip()+'\tFalse negative\n')
		#True negative
		sum_of_all = (int(Genome_completeness) + int(Resultat_du_test_rapide) + int(Reultat_PCR_gargarisme) + int(Resultat_PCR_nasal) + int(Postive_RT_PCR) + MGC_qPCR_results2)
		u = 1 if (MGC_qPCR_results == 'Negative' and sum_of_all == 0) else 0
		v = 1 if Genome_completeness1 == '0' and sum_of_all == 0 else 0
		w = 1 if (Resultat_du_test_rapide1 == '0' and sum_of_all == 0) else 0
		x = 1 if (Reultat_PCR_gargarisme1 == '0' and sum_of_all == 0) else 0
		y = 1 if (Resultat_PCR_nasal1 == 0 and sum_of_all == 0) else 0
		z = 1 if (Postive_RT_PCR1 == 0 and sum_of_all == 0) else 0
		if u == 1 or v == 1 or w == 1 or x == 1 or y == 1  or z == 1:
			if lamp_result == 'Positive':
				newfile.write(line.strip()+'\tFalse positive\n')
			elif lamp_result == 'Negative':
				newfile.write(line.strip()+'\tTrue negative\n')
	newfile.close()
	
def corinne_file_man3(corinne_file):
	newfile = open(corinne_file+'_ed2','w')
	for line in open(corinne_file, 'r'):
		gtruth = line.split()[0]
		lamp = line.split()[1]
		score = 'NA'
		if gtruth == 'Negative':
			if lamp == gtruth:
				score = 'TN'
			elif lamp == 'Positive':
				score = 'FP'
			else:
				score = 'NA'
		elif gtruth == 'Positive':
			if lamp == gtruth:
				score = 'TP'
			elif lamp == 'Negative':
				score = 'FN'
			else:
				score = 'NA'
		else:
			score = 'NA'
		newfile.write(line.strip()+'\t'+score+'\n')
	newfile.close()
	
def return_genome_coverage(samtools_depth):
	newfile = open(samtools_depth+'_ed2','w')
	
	dic = {}
	filey = open(samtools_depth, 'r')
	filex = filey.readlines()
	last_pos = filex[-1].split()[1]
	ref = filex[-1].split()[0]		
	print(last_pos)
	
	for a in range(len(filex)):
		linex = filex[a].split('\t')
		ref = linex[0]
		pos = linex[1]
		cov = linex[2]
		dic.setdefault(ref+'_'+pos,int(cov))
		for b in range(1,int(pos)+1):
			if ref+'_'+str(b) not in dic.keys():
				dic.setdefault(ref+'_'+str(b),0)
	for c in range(1,29903):
		if ref+'_'+str(c) not in dic.keys():
			dic.setdefault(ref+'_'+str(c),0)
	for ref_pos,cov in dic.items():
		newfile.write(ref_pos.split('_')[0] + '\t' + ref_pos.split('_')[1] + '\t' + str(cov) + '\n')
	newfile.close()
	
def return_genome_coverage2(samtools_depth):
	newfile = open(samtools_depth+'_ed2','w')	
	dic = {}
	filey = open(samtools_depth, 'r')
	filex = filey.readlines()
	for i in range(0,29903):
		pos1 = i+1
		cov = 0
		for linex in filex:
			ref = linex.split('\t')[0]
			pos = linex.split('\t')[1]			
			if int(pos) == pos1:
				cov = linex.split('\t')[2]
		dic.setdefault(ref+'_'+str(pos1),int(cov))
		
	for ref_pos,cov in dic.items():
		newfile.write(ref_pos.split('_')[0] + '\t' + ref_pos.split('_')[1] + '\t' + str(cov) + '\n')
	newfile.close()
	
def return_genome_mean_cov(rootdir):
	import statistics as stats
	dic = {}
	samtools_depth_dir = '/home/banthony/scratch/analysis/covid19/B004_11_2/plot_amplicons/depth_files'
	for dirName, subdirList, fileList in os.walk(rootdir):
		if len(subdirList) == 0 and len( [x for x in fileList if x.endswith('bam')] ) > 0:
			bams = [x for x in fileList if x.endswith('bam') ]
			for bam_file in bams:
				protocol = bam_file.split('-')[0]
				variant = bam_file.split('-')[1]
				reads = bam_file.split('-')[2].replace('.bam','')
				samtools_depth_file = bam_file.replace('bam','samtools_depth')
				samtools_depth_file = samtools_depth_dir+'/'+samtools_depth_file
				#os.system("samtools depth -aa %s > %s" %(dirName+'/'+bam_file, samtools_depth_file) )
				os.system("samtools depth %s > %s" %(dirName+'/'+bam_file, samtools_depth_file) )
				cov_list = []
				for line in open(samtools_depth_file, 'r'):
					cov = line.strip().split('\t')[2]
					if int(cov) >= 20:
						cov_list += [int(cov)]
					else:
						cov_list += [0]
				stdev = str(round(stats.stdev(cov_list), 3))
				mean = str(round(stats.mean(cov_list), 3))
				try:
					cov = str( round((100*float(stdev)/float(mean)),3))
				except:
					cov = '0'
				dic.setdefault( '\t'.join([protocol,variant,reads]), cov+'\t'+mean+'\t'+stdev+'\n' )
	
	newfile = open(rootdir+'/combined_depths_stats3','w')
	newfile.write('\t'.join(['protocol','variant','reads','CoV','mean','stdev\n']))
	for x,y in dic.items():
		newfile.write('\t'.join([x,y]))
	newfile.close()

	
def look_up_given_reads_fasta(file,fasta): #Use this for speed
	#This script takes a fasta file and a file with the names of reads u want to take from the fasta. It will create a dictionary for the fasta seqs and search read names against the fasta name. 
	#If a name in the file is found in the fasta as well, the corresponding sequence and the name are returned. Pay attention to the \n and > and fasta format, it has to be single line for seqs
	
	print("\nCreating dictionary of the fasta...")
	#lopenfasta = fasta_to_dict(fasta)
	lopenfasta = fastq_to_dict(fasta)
	openfile=open(file,'r')
	lopenfile = openfile.readlines()
	lopenfile = list(set(lopenfile))
	l = len(lopenfile)
	#path = '/home/banthony/scratch/reads/nanopore/covid19/B004/B004_twist1/mapped_reads'
	outa=open(path+'/'+fasta.strip('/').split('/')[-1],'w')
	
	#Check
	print([lopenfile[0]])
	print(list(lopenfasta.keys())[0])
	total_reads = len(lopenfasta.keys())
	a=b=c=d=e=f=n = 0
	print("\nQuerrying fasta dictionary for reads...\n")
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for read in lopenfile:
		readnamex = read.split()[0].strip() #+'\n'			
		if readnamex in lopenfasta.keys(): #'>'+read.strip().split()[0]+'\n' in lopenfasta.keys():
			outa.write(lopenfasta[readnamex])
			a += 1
		n+=1
		bar.update(n)
	bar.finish()	
	outa.close()
	print('\nDiscovered ' + str(a) + ' sequences' + ' which is %s percent of all the names u gave' %(str(round(100*a/(l-2)))))
	
def filter_fastq_rdlen(fasta):
	outa=open(fasta+'_ed','w')
	lopenfasta = fastq_to_dict(fasta)
	for x,y in lopenfasta.items():
		if len(y.split('\n')[1]) > 199:
			outa.write(y)
	outa.close()
	
def select_aligned_reads(rootdir):
	global path
	dirs = ['artic_v3','artic_v4','artic_v4.1','entebbe','midnight','snap','1','1e2','1e3','1e4','1e5','1e6','10']
	ref = '/home/banthony/projects/rrg-ioannisr/banthony/analysis/covid19/genome/NC_045512.2.fasta'
	mapped_reads = path = '/home/banthony/scratch/reads/nanopore/covid19/B004/B004_twist1/mapped_reads/'
	for dirName, subdirList, fileList in os.walk(rootdir):
		if len(subdirList) == 0 and dirName.strip().split('/')[-2] in dirs:			
			fastq = [x for x in fileList if x.endswith('.fastq') and not 'pchop1' in x]
			for i in fastq:
				print(dirName+'/'+i)
				#align reads and return aligned reads
				destin1 = dirName.replace('B004_twist1x','sra')
				os.makedirs(destin1,exist_ok=True)
				mapped_readsx = mapped_reads+i.replace('.fastq','.mapped_reads')
				os.system(' minimap2 -a --sam-hit-only -t 23 %s %s | cut -f 1 > %s' %(ref, dirName+'/'+i, mapped_readsx))
				
				openfilea=open(mapped_readsx,'r')
				lopenfilea = openfilea.readlines()
				if len(lopenfilea) < 5:
					continue
				else:				
					#purse fastq and return mapped reads
					look_up_given_reads_fasta(mapped_readsx,dirName+'/'+i)
					os.system('mv %s %s' %(mapped_reads+'/'+i,destin1))
					
def select_aligned_reads2(rootdir):
	global path
	ref = '/home/banthony/projects/rrg-ioannisr/banthony/analysis/covid19/genome/NC_045512.2.fasta'
	mapped_reads = path = '/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/mapped_reads/'
	for dirName, subdirList, fileList in os.walk(rootdir):
		if len(subdirList) == 0:
			if 'entebbe' in dirName or 'midnight' in dirName or 'snap' in dirName or 'artic_v3' in dirName or 'artic_v4.1' in dirName:			
				#fastq = [x for x in fileList if x.endswith('.fastq') and not 'pchop1' in x]
				fastq = [x for x in fileList if x.endswith('.fastq')]
				for i in fastq:
					print(dirName+'/'+i)
					#align reads and return aligned reads
					destin1 = dirName.replace('fastq_pass','sra')
					os.makedirs(destin1,exist_ok=True)
					mapped_readsx = mapped_reads+i.replace('.fastq','.mapped_reads')
					#os.system(' minimap2 -a --sam-hit-only -t 23 %s %s | cut -f 1 > %s' %(ref, dirName+'/'+i, mapped_readsx))

					openfilea=open(mapped_readsx,'r')
					lopenfilea = openfilea.readlines()
					if len(lopenfilea) < 5:
						continue
					else:				
						#purse fastq and return mapped reads
						look_up_given_reads_fasta(mapped_readsx,dirName+'/'+i)
						os.system('mv %s %s' %(mapped_reads+'/'+i,destin1))
						
def select_aligned_reads3(rootdir):
	global path
	dirs = ['artic_v3','artic_v4','artic_v4.1','entebbe','midnight','snap','1','1e2','1e3','1e4','1e5','1e6','10']
	ref = '/home/banthony/projects/rrg-ioannisr/banthony/analysis/covid19/genome/NC_045512.2.fasta'
	mapped_reads = path = '/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_3/mapped_reads/'
	for dirName, subdirList, fileList in os.walk(rootdir):
		fastq = [x for x in fileList if x.endswith('.fastq') and not 'pchop1' in x and 'barcode' in x]
		if len(fastq) > 0 and 'R2C2' not in dirName:
			if 'OM' in dirName or 'IND' in dirName:
				for i in fastq:
					print(dirName+'/'+i)
					#align reads and return aligned reads
					destin1 = dirName.replace('fastq_pass','sra')
					os.makedirs(destin1,exist_ok=True)
					mapped_readsx = mapped_reads+i.replace('.fastq','.mapped_reads')
					os.system(' minimap2 -a --sam-hit-only -t 23 %s %s | cut -f 1 > %s' %(ref, dirName+'/'+i, mapped_readsx))

					openfilea=open(mapped_readsx,'r')
					lopenfilea = openfilea.readlines()
					if len(lopenfilea) < 5:
						continue
					else:				
						#purse fastq and return mapped reads
						#print([mapped_readsx,dirName+'/'+i])
						look_up_given_reads_fasta(mapped_readsx,dirName+'/'+i)
						os.system('mv %s %s' %(mapped_reads+'/'+i,destin1))
						
def select_unligned_reads4(rootdir):
	from shlex import quote
	ref = '/home/banthony/scratch/analysis/nanopore/capitata/ncbi_genome/genome_ercc/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_plus_EGII_novelgenes.fasta.minimap2.idx'
	mapped_reads = rootdir+'/unmapped/'
	os.makedirs(mapped_reads,exist_ok=True)
	mapped_reads = rootdir+'/unmapped/unmapped.fasta'
	mapped_reads2 = rootdir+'/unmapped/unmapped2.fasta'
	all_reads = []
	newfile = open(mapped_reads, 'w')
	for dirName, subdirList, fileList in os.walk(rootdir):
		if len(subdirList) == 0 and 'barcode' not in dirName and 'unclassified' not in dirName:
			#fastq = [x for x in fileList if x.endswith('all.porechop.pychop.cutadapt.fq') ]
			fastq = [x for x in fileList if x.endswith('.fastq') or x.endswith('.fasta') ]
			for i in fastq:
				print(dirName+'/'+i); all_reads += [dirName+'/'+i]
				os.system(" minimap2 -ax splice -t 23 --secondary=no %s %s | awk ' $2==4 || $3 == \"*\" {print(\">\"$1,$10)}' > %s" %(ref, dirName+'/'+i, mapped_reads2))				
				for linex in open(mapped_reads2, 'r'):
					if len(linex.split()[1]) > 100:
						newfile.write('\n'.join(linex.split())+'\n')
				#os.system(" minimap2 -a -t 23 %s %s | awk '$2==4 {print(\">\"$1\"\n\"$10)}' >> %s" %(ref, dirName+'/'+i, mapped_reads))
				#cmd = '$2==4 || $3 == "*" {print(">"$1"\n"$10)}'
				#os.system("minimap2 -a -t 23 %s %s | awk %s >> %s" %(ref, dirName+'/'+i, quote(cmd), mapped_reads))
				#os.system(" minimap2 -a -t 23 %s %s | awk '$2==4 || $3 == \"*\" {print(\">\"$1\"\n\"$10)}' >> %s" %(ref, dirName+'/'+i, mapped_reads))
	newfile.close()
	ref = '/home/banthony/scratch/analysis/nanopore/capitata/ncbi_genome/CAJHJT01/CAJHJT01.1.fsa_nt.minimap2.idx'
	wkdir = rootdir+'/unmapped/'
	output = 'unmapped'
	os.system(' cd %s && minimap2 -ax splice -t 23 --secondary=no %s %s > %s' %(wkdir, ref, mapped_reads, mapped_reads.replace('fasta','sam')))
	os.system(' cd %s && samtools view --threads 23 -b -o %s %s' %(wkdir, output+'.bam', output+'.sam'))
	os.system(' cd %s && samtools sort --threads 23 -o %s %s' %(wkdir, output+'.sort.bam', output+'.bam'))
	os.system(' cd %s && samtools index %s' %(wkdir, output+'.sort.bam'))
	os.system(' cd %s && rm %s' %(wkdir, output+'.bam'))
	print('\n'.join(all_reads))
	
def collect_align_fastq(rootdir, output):
	ref = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/unmapped/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_plus_EGII_novelgenes2.fasta'
	annot1 = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_FINAL_flair.gtf'
	annot2 = '/home/banthony/scratch/analysis/nanopore/capitata/ncbi_genome/annotation/ERCC92_GCF_000347755.3_Ccap_2.1_genomic.gtf'
	'''
	for dirName, subdirList, fileList in os.walk(rootdir):
		if len(subdirList) == 0 and 'barcode' not in dirName and 'unclassified' not in dirName:
			fastq = [x for x in fileList if x.endswith('all.porechop.pychop.cutadapt.fq') ]
			for i in fastq:
				print(dirName+'/'+i)
				os.system('cat %s >> %s' %(dirName+'/'+i, rootdir+'/'+'combined.all.porechop.pychop.cutadapt.fq'))
	'''
	mapped_reads = rootdir+'/'+'combined.all.porechop.pychop.cutadapt.fq'
	wkdir = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/'
	#output = 
	'''
	os.system(' cd %s && minimap2 -ax splice -t 23 --secondary=no %s %s > %s' %(wkdir, ref, mapped_reads, output+'.sam'))
	os.system(' cd %s && samtools view --threads 23 -b -o %s %s' %(wkdir, output+'.bam', output+'.sam'))
	os.system(' cd %s && samtools sort --threads 23 -o %s %s' %(wkdir, output+'.sort.bam', output+'.bam'))
	os.system(' cd %s && samtools index %s' %(wkdir, output+'.sort.bam'))
	os.system(' cd %s && rm %s' %(wkdir, output+'.bam'))'''
	os.system(' cd %s && stringtie -j 5 -c 10 -s 10 -p 23 -G %s -l %s -A %s -L -o %s %s' %(wkdir, annot2, output+'_annot2.stringtie.transcripts', output+'_annot2.stringtie.quantification', output+'_annot2.stringtie.gtf', output+'.sort.bam'))
	
def run_stringtie(rootdir):
	ref = '/home/banthony/scratch/analysis/nanopore/capitata/ncbi_genome/genome_ercc/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_plus_EGII_novelgenes.fasta.minimap2.idx'
	annot1='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_flair_tofu_stringtie.sorted.gtf'

	base_dir = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/'
	for dirName, subdirList, fileList in os.walk(rootdir):
		if len(subdirList) == 0 and 'barcode' not in dirName and 'unclassified' not in dirName and 'unmapped' not in dirName:
			fastq = [x for x in fileList if x.endswith('all.porechop.pychop.cutadapt.fq') ]
			for i in fastq:
				print(dirName+'/'+i)
				same_id = output = i.replace('.all.porechop.pychop.cutadapt.fq','')
				same_id2 = same_id[0:2]+'H' if len(same_id) == 3 else same_id[0:3]+'H'
				wk_dir = base_dir+same_id2+'/'+same_id+'/'
				os.makedirs(wk_dir,exist_ok=True)
				os.system(" cd %s && rm *stringtie*" %(wk_dir))
				#os.system(" cd %s && minimap2 -ax splice -t 23 --secondary=no %s %s > %s" %(wk_dir,ref, dirName+'/'+i, same_id+'.sam'))
				#os.system(' cd %s && samtools view --threads 23 -h -b -q 0 -o %s %s' %(wk_dir, output+'.bam', output+'.sam'))
				#os.system(' cd %s && samtools sort --threads 23 -o %s %s' %(wk_dir, output+'.sort.bam', output+'.bam'))
				#os.system(' cd %s && samtools index %s' %(wk_dir, output+'.sort.bam'))
				#os.system(' cd %s && rm %s' %(wk_dir, output+'.bam'))
				os.system(' cd %s && stringtie --fr -e -p 23 -G %s -A %s.stringtie.tpm -L -o %s.stringtie.gtf %s.sort.bam' %(wk_dir, annot1, same_id, same_id, same_id))
				#os.system(' mv %s %s' %(wk_dir+same_id+'.stringtie.tpm','/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie_expn/all/'+same_id+'.stringtie.tpm'))
				os.system(' cp %s %s' %(wk_dir+same_id+'.stringtie.tpm','/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie_expn/all/'))
				
def return_read_names(rootdir):
	listdir = []
	for dirName, subdirList, fileList in os.walk(rootdir):
		if len(subdirList) == 0:
			for filex in fileList:
				listdir += ['\t'.join([dirName,'\t'.join(dirName.strip('/').split('/')), filex+'\n'])]
	outa=open(rootdir+'/files','w')
	outa.write(''.join(listdir))
	outa.close()
	
def combine_fastq_files(rootdir, details):
	det = {}
	newfile = open(rootdir+'/Cell_culture_experiment.fastq','w')
	for linex in open(details, 'r'):
		det.setdefault(linex.split()[0], linex.split()[1])
	for dirName, subdirList, fileList in os.walk(rootdir):
		if len(subdirList) == 0 and len([x for x in fileList if x.endswith('.fastq')]) > 0:
			print("Processing " + dirName +', '+str(len(fileList))+' files')
			for filex in fileList:
				if filex in det.keys():
					lopenfasta = fastq_to_dict(dirName+'/'+filex)
					for x,y in lopenfasta.items():
						y = y.lstrip('@')
						newfile.write('@'+det[filex]+':'+y)
	newfile.close()
	
def rename_fastq_readsx(files):
	for filex in files:
		newfile = open(filex+'_ed','w')
		for line in open(filex, 'r'):
			if line.startswith('@'):
				line = line.replace('ARTIC_v', 'ARTIC-v')
				line = line.replace('Twist_Control_1','Twist-Control-1')
				newfile.write(line)
			else:
				newfile.write(line)
	newfile.close()

def combine_stringtie_quantification(rootdir):
	all_dic = {}
	x = 0; print("Creating dictionary...\n")
	bar = progressbar.ProgressBar(maxval=len(os.listdir(rootdir)), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()	
	for filex in os.listdir(rootdir):
		sample_name = filex.split('.')[0]
		all_dic.setdefault(sample_name,{})
		for linex in open(rootdir+'/'+filex, 'r'):
			if not linex.startswith('Gene'):
				gene_id = linex.split()[0]
				all_dic[sample_name].setdefault(gene_id,linex)
		x += 1
		bar.update(x)
	bar.finish()		
	genes = sorted(list(all_dic[list(all_dic.keys())[0]].keys()))
	samples = sorted(list(all_dic.keys()))
	
	combined_results_file = rootdir.rstrip('/').rstrip(rootdir.strip('/').split('/')[-1])+'/C0_15H.combined.stringtie.tpm'; os.system('rm %s' %combined_results_file)
	newfile = open(combined_results_file,'w')
	newfile.write('\t'.join([ 'Gene_ID','Gene_name','Reference','Strand','Start','End' ]) + '\t') #Remove 'Gene_Name', since its mostly empty
	newfile.write('\t'.join(samples) + '\n')
	
	x = 0; print("\nPursing dictionary...")	
	bar = progressbar.ProgressBar(maxval=len(genes), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for genex in genes:
		Gene_Name = all_dic[samples[0]][genex].split()[1]
		gene_info = '\t'.join(all_dic[samples[0]][genex].split()[0:6])
		newfile.write( gene_info + '\t') #newfile.write( gene_info.replace(Gene_Name,'').replace('\t\t','\t') + '\t')
		for samplex in samples:
			#print([samplex, genex]) #print([samplex, genex, all_dic[samplex][genex]])
			newfile.write(  all_dic[samplex][genex].split()[8]  + '\t')
		newfile.write('\n')
		x += 1
		bar.update(x)
	bar.finish()
	newfile.close()
	'''
	all_samples = []
	for filex in os.listdir(rootdir):
		sample_name = filex.split('.')[0]
		all_samples += [sample_name]
		os.system('cd %s && sort -k1,1 %s | cut -f 9 | paste %s - > %s' %(rootdir, rootdir+'/'+filex, 'C0_15H.combined.stringtie.tpm'))
		#os.system('paste ')
	'''	
def combine_rsem_to_stringtie(rsem, stringtie):
	rsem_dic = {}
	a=0; b=0
	for line in open(rsem, 'r'):
		gene = line.split()[0]
		if gene not in rsem_dic.keys():
			rsem_dic.setdefault(gene, line)
		else:
			print(gene)
	newfile = open(stringtie+'.rsem.tpm','w')
	for line in open(stringtie, 'r'):		
		if a == 0:
			newfile.write(line.strip()+'\trsem_tpm\tmean\n')			
		else:
			#print(line.split()[:5])
			if 'ERCC' in line and 'ERCC-' not in line:
				line = line.replace('ERCC','ERCC-')
			gene = line.split()[0]
			gene_mean_ls = line.strip().split()[6:]
			for i in gene_mean_ls:
				b += float(i)
			gene_mean = str(round((b/len(gene_mean_ls)), 3))
			b = 0
			try:				
				rsem_tpm = rsem_dic[gene].split()[5]				
			except:
				rsem_tpm = '0'
			newfile.write(line.strip()+'\t'+'\t'.join([rsem_tpm, gene_mean]) + '\n')
		a += 1
	newfile.close()
	
def multiply_reads(fasta):
	lopenfasta = fasta_to_dict(fasta)	
	newfile = open(fasta+'2','w')
	for name,seq in lopenfasta.items():
		a = 1
		for i in range(0,20,1):
			newfile.write(name.strip()+'_'+str(a)+'\n'+seq)
			a += 1
	newfile.close()
	
def flair_read_matrix(rootdir, sampleid_sex, batch):
	sample_dic = {}
	for line in open(sampleid_sex, 'r'):
		sample_dic.setdefault(line.split()[1],line)
	
	newfile = open('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/flair/'+'reads_manifestC.txt','w')

	for dirName, subdirList, fileList in os.walk(rootdir):
		if len(subdirList) == 0 and 'barcode' not in dirName and 'unclassified' not in dirName and 'unmapped' not in dirName and len([x for x in fileList if x.endswith('all.porechop.pychop.cutadapt.fq') ])> 0:
			fastq = [x for x in fileList if x.endswith('all.porechop.pychop.cutadapt.fq') ]
			for i in fastq:
				sample_dets = sample_dic[i.split('.')[0]].split()
				sample_id = sample_dets[1]
				sex = sample_dets[2]
				sample_id = sample_id.replace(sample_id[-1],'H'+sample_id[-1])+'.'+sex
				condition = 'conditionA' if sex == 'X' or sex == 'F' else 'conditionB'
				location = dirName+'/'+i
				newfile.write('\t'.join([sample_id,condition, batch,location])+'\n')
				
def add_average_tpm(classn, flair, TPM_cutoff):
	file_ext = '.'+classn.split('.')[-1] if classn.split('.')[-1] != "."  else classn[-3:-1]
	flair_dic = {}
	for line in open(flair, 'r'):
		flair_dic.setdefault(line.split('\t')[0], line.strip().split('\t')[1:])
	newfile = open(classn.replace(file_ext, '_ed'+file_ext), 'w'); b=0
	for line in open(classn, 'r'):
		if line.startswith('isoform'):
			newfile.write(line.strip()+'\tperc_samples_with_TPM_above'+str(TPM_cutoff)+'\n' )
		else:
			a=0
			isoform = line.split('\t')[0]
			try:
				tpm_ls = flair_dic[isoform]
			except:
				tpm_ls = [0]
				b += 1
			for i in range(len(tpm_ls)):			
				if float(tpm_ls[i]) >= TPM_cutoff:
					a += 1
			rep = str(round(100*a/len(tpm_ls),3))
			if len(line.split('\t')) != 551:
				newfile.write(line.strip()+'\tND\t'+rep+'\n')
			else:
				newfile.write(line.strip()+'\t'+rep+'\n')
	newfile.close(); print(b)

def get_longest_orf(faa, bed):
	longest_orf = []
	for line in open(bed, 'r'):
		transcript = line.split()[0]
		pattern = re.compile(r'_ORF.(\d+);')
		orf = '>' + transcript + '_ORF.'+ pattern.search(line).group(1)
		longest_orf += [orf]
	reads_fasta = check_fasta.check_fasta_fmt(faa)
	newfile = open(faa.replace('faa', 'longest.faa'), 'w')
	length = int(os.popen('wc -l %s' %(reads_fasta)).read().split()[0])
	file1 = open(reads_fasta,'r'); k=0
	while k < length:
		read_name =file1.readline()
		seq =file1.readline()
		if read_name.split()[0] in longest_orf:
			newfile.write(read_name+seq)
		k+= 2
	newfile.close()
	
def reorder_columns(file1, correct_order): #This script correctly orders the colums of gene expression files generated from Stringtie or Flair.
	#In case you have a file with the correct way to order the columns
	for line in open(correct_order, 'r'):
		correct_order_ls = [k.strip('\'') for k in line.strip('[]').split(', ')]
	#Just put all lines of file1 in a dictionary
	items_dict = {}
	for line in open(file1, 'r'):
		items_dict.setdefault(line.split()[0], line)
	#This is the heavy code. It basically stores all column entries in one dictionary list. So, every column is assigned a key which is the column
	#label and value which is a list of the column label and that column's entries.
	items_dict2 = {}; orderx_dict = {}; x = 0
	for line in open(file1, 'r'):
		line_fields = line.strip().split()
		if x == 0:			
			for field_indx in range(len(line_fields)):
				items_dict2.setdefault(line_fields[field_indx], [line_fields[field_indx]])
				orderx_dict.setdefault(field_indx, line_fields[field_indx])
		else:
			for field_indx in range(len(line_fields)):
				items_dict2[orderx_dict[field_indx]].append(line_fields[field_indx])
		x+=1
	#Just get the header
	try:
		head = items_dict['Gene_ID'].strip().split()
	except:
		head = items_dict['ids'].strip().split()
	
	#Create the correct order of columns
	timepoint_pat = re.compile(r'C(\d+)\w')	
	order = []
	for i in head:
		if not i.startswith('C'):
			order += [i]
	for i in range(0,16,1):
		for sample in [x for x in head if x.startswith('C')]:
			timepoint = timepoint_pat.search(sample).group(1)
			if str(i) == timepoint and sample != 'C0_15H':
				order += [sample]
	#print(['order'] + order)		
	#Write the new ordered file
	newfile = open(file1+'_ed', 'w')
	for i in range(len(items_dict2[list(items_dict2.keys())[0]])):
		to_write = []
		for column in order:
			to_write += [items_dict2[column][i]]
		newfile.write('\t'.join(to_write)+'\n')
	newfile.close()
	
def change_stringtie_colnames(stringtie, correctnames):
	for line in open(correctnames, 'r'):
		line = line.strip().strip('[]')
		correct_order_ls = [k.strip('\'') for k in line.strip('[]').split(', ')]
	newfile = open(stringtie+'_ed', 'w'); k=0
	for line in open(stringtie, 'r'):
		if k == 0:
			new_names = []
			line_items = line.strip().split()
			for i in line_items:
				if i.startswith('C'):
					sample_pat = re.compile(r'C\d+(\w)')
					sample = sample_pat.search(i).group(1)
					new_sample = i.rstrip(sample)+ 'H'+'.'+sample+'.'
					for j in correct_order_ls:
						if new_sample in j:
							new_names += [j.strip('[]').strip('\'')]
				else:
					new_names += [i]
			newfile.write('\t'.join(new_names) + '\n'); k+=1
		else:
			newfile.write(line)
	newfile.close()
	
def calc_gene_expn(file1):
	items_dict2 = {}; orderx_dict = {}; xx=yy = 0
	for line in open(file1, 'r'):
		line_fields = line.strip().split()
		if xx == 0:
			line_fields2 = ['Gene_ID.X']
			for idx in line_fields:
				sample_id = idx.split('.')[0] + '.' + idx.split('.')[2]
				line_fields2 += [sample_id]				
			for field_indx in range(len(line_fields2)):
				sample_id = line_fields2[field_indx]
				if sample_id not in items_dict2.keys():
					items_dict2.setdefault(sample_id, {})
				orderx_dict.setdefault(field_indx, sample_id)
		else:
			gene_id = line_fields[0]
			for x,y in items_dict2.items():
				y.setdefault(gene_id, [])
			for field_indx in range(len(line_fields)):
				val = line_fields[field_indx]
				if field_indx == 0:
					continue #items_dict2[orderx_dict[field_indx]][gene_id].append(val)
				else:
					try:
						val2 = float(val)
						items_dict2[orderx_dict[field_indx]][gene_id].append(val2)
						if val2 <= 5000:
							yy += 1
					except:
						kk = 0 #print(line.split())
						#sys.exit(gene_id + ' is not a float ' + str(val))
		xx +=1
	#print(items_dict2[ 'Gene_ID.X'])
	samples = sorted(list(items_dict2.keys()))
	samples = [x for x in samples if not x.startswith('Gene_ID')]
	genes = sorted(list(items_dict2['Gene_ID.X'].keys()))
	genes = [x for x in genes if not x.startswith('ERCC')]
	newfile = open(file1+'_ed', 'w')
	newfile.write( 'Gene_ID'+'\t'+'\t'.join(samples) + '\n' )
	for gene in genes:
		averages = []
		for sample in samples:
			exp_ls = items_dict2[sample][gene]
			#lenx = len([m for m in exp_ls if m != 0]) if len([m for m in exp_ls if m != 0]) != 0 else 1
			lenx = 1 if sum(exp_ls) == 0 else len([m for m in exp_ls if m != 0])
			average = round(sum(exp_ls)/lenx, 2)
			averages += [str(average)] #averages += [str(average).split('.')[0]]
		newfile.write( gene+'\t'+'\t'.join(averages) + '\n')
	newfile.close(); print([yy, 'Done....'])
	
dr = '/home/banthony/scratch/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/SQANTI_tofu_min2_polyAtrim/novel_genes/minimap2_find_wrong_novel_genes'
dr = '/home/banthony/scratch/analysis/nanopore/olivefly/ncbi_genome/full_annotation/proteins/functional_annotation'
dr = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/tofu'
dr = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/tofu/reflair/sqanti/side_job'

#manip_gbk(dr + '/' + 'protein.gbk')
#manip_gbk(dr + '/' + '12')
#remove_bad_samlines2(dr+'/'+'picard.stdout', dr+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs_Bo_E_6H_C010_07_pass_edited.strand_no_strand_choped_sl.subsample_2M_out_1_minimap2.sam')
#dist_reads(dr+'/'+'Bo.E.Heads.pass.corrected.stranded.choped.modified.cutadapt.fasta')
#get_metrichor_minus_nanonet('/gs/project/wst-164-ab/anthony/Nanopore_data/ecoli/ecoliTest/C002_05_6_50kb/metrichor/9.txt')
#find_gtf_geneID(dr+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_gffread_sorted_edited_over37.gtf2')
#get_ercc_mapping(dr+'/'+'ercc_read_names', dr+'/'+'ERCC_mapping_filtered_reads.txt_2_readlengths_GC.txt', dr+'/'+'ERCC92_sl_readlengths_GC.txt', dr+'/'+'ERCC_mapping_filtered_reads2.txt')
#smthing(dr+'/'+'ERCC_ratios',dr+'/'+'ERCC92_sl_readlengths_GC.txt')
#ERCC92_sl_readlengths_GC.txt
#ERCC_mapping_filtered_reads.txt_2_readlengths_GC.txt
#ERCC_mapping_filtered_reads2.txt
#add_rdlen_readname(dr+'/'+'drosophila_melanogaster_uniprot-proteome_UP000000803.fasta')
#add_rdlen(dr+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds.rna_sl_2_Bo_E_1H_C010_10_pass_edited.trimmed_stranded_edited_gmap_with_alignments.txt')
#change_read_names(dr+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs.fasta',dr+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds.gff3.gpd')
#manip_gff3(dr+'/'+'GCF_000347755.3_Ccap_2.1_genomic.gff') #z.get.transcripts')
#manip_gff3(dr+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds.gff3')
#gff3_to_gtf(dr+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gff3_1')
#gff3_to_gtf_ftr(dr+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds.gff3')
#find_pattern(dr+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gff3_1')
#correct_files2(dr+'/'+'Bo_E_all_pass_edited.correctedReads.no_adapters_choped.fasta')
#correct_files2(dr+'/'+'Bo_E_all_pass_edited.correctedReads.single_adapter_choped.fasta')
#correct_files2(dr+'/'+'Bo_E_5H_pass_edited_canu_corrected.choped.fasta')
#correct_files3(dr+'/'+'Bo_E_all_pass_edited.correctedReads.trimmed_stranded_choped.fasta')
#correct_files(dr+'/'+'Bo_E_1H_C010_10_pass_edited.no_adapters.fasta')
#correct_files(dr+'/'+'Bo_E_1H_C010_10_pass_edited.single_adapter.fasta')
#correct_files(dr+'/'+'Bo_E_2H_C010_09_pass_edited.no_adapters.fasta')
#correct_files(dr+'/'+'Bo_E_2H_C010_09_pass_edited.single_adapter.fasta')

#correct_files(dr+'/'+'Bo_E_3H_C010_08_pass_edited.no_adapters.fasta')
#correct_files(dr+'/'+'Bo_E_3H_C010_08_pass_edited.single_adapter.fasta')

#correct_files(dr+'/'+'Bo_E_4H_C010_06_pass_edited.no_adapters.fasta')
#correct_files(dr+'/'+'Bo_E_4H_C010_06_pass_edited.single_adapter.fasta')

#correct_files(dr+'/'+'Bo_E_5H_pass_edited.no_adapters.fasta')
#correct_files(dr+'/'+'Bo_E_5H_pass_edited.single_adapter.fasta')

#correct_files(dr+'/'+'Bo_E_all_pass_edited.correctedReads.no_adapters.fasta')
#correct_files(dr+'/'+'Bo_E_all_pass_edited.correctedReads.single_adapter.fasta')
#dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_4H/C010_06_191017/mandalorion'
#get_reverse(dr + '/' + '2D_trimmed_l_filtered.fasta', 'fasta')
#dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_5H/combined/mandalorion'
#get_reverse(dr + '/' + '2D_trimmed_l_filtered.fasta', 'fasta')
#dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_6H/C010_07_4_041117/mandalorion'
#get_reverse(dr + '/' + '2D_trimmed_l_filtered.fasta', 'fasta')
#manip_samfile(dr + '/' + 'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs_10X_only_gene_refs_Bo-E-6H_illumina.noERCC.sam')
#manip_samfile(dr + '/' + 'header')
#manip_samfile(dr + '/' + '123.sam')
#manip_samfile(dr + '/' + 'taylorFusions_C012_01_mb704.passFail_minimap.sam')
#manip_samfile(dr + '/' + 'taylorFusions_C012_01_mb2018.passFail_minimap.sam')
#manip_samfile(dr + '/' + 'longest_transcript_6H.male.female_minimap2.sam')

#find_uniprot2(dr+'/'+'uniprot_trembl_headers',dr+'/'+'no_uniprot_sprot_annnotations_exact_match')
#find_uniprot2(dr+'/'+'uniprot_trembl_headers',dr+'/'+'no_uniprot_sprot_annnotations_exact_match_after_removing-like') #GO2_sorted') #zv2') #GO2_sorted') #
#find_uniprot2(dr+'/'+'uniprot_trembl_annnotations_exact_match_headers',dr+'/'+'no_uniprot_sprot_annnotations_exact_match_after_removing-like')
#find_uniprot2(dr+'/'+'uniprot_trembl_headers',dr+'/'+'no_uniprot_sprot_annnotations_exact_match2')
#find_common_lines2(dr+'/'+'no_uniprot_sprot_annnotations_exact_match2', dr+'/'+'uniprot_trembl_annnotations_exact_match_headers_formated_aa_grep_unique')
#find_common_lines2(dr+'/'+'representative_seqs_for_10X_genes', dr+'/'+'best.sorted_resorted.bed')
#dist_reads2(dr+'/'+'contig131784_reverse_Bo_EnHeads_pass_both_adapters_stranded.choped.manda_trim_minimap2_ed.sam')
#dist_reads2(dr+'/'+'Bo_E_all_pass_edited_intronicreads.Bo.E.Heads')
#dist_reads3(dr+'/'+'1/ncbi_genome_10X_only_gene_refs_Bo.E.Heads.pass.corrected.stranded.choped.modified.cutadapt_gmap.sorted.noERCC.sam')
#change_read_names_to_pacbio(dr + '/' + 'C1_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_15/C1H/pass/cutadapt/C1_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_14/C2H/pass/cutadapt/C2_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_13/C3H/pass/cutadapt/C3_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_12/C4H/pass/cutadapt/C4_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_11/C5H/pass/cutadapt/C5_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_10/C6H/pass/cutadapt/C6_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_9/C7H/pass/cutadapt/C7_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_8/C8H/pass/cutadapt/C8_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_7/C9H/pass/cutadapt/C9_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_6/C10H/pass/cutadapt/C10_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_5/C11H/pass/cutadapt/C11_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_4/C12H/pass/cutadapt/C12_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_3/C13H/pass/cutadapt/C13_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/reads/nanopore/capitata/C011_9_2/C14H/pass/cutadapt/C14_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
#change_read_names_to_pacbio('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/tofu/reflair/sqanti/transcript_ids_for_flair_isoforms.fasta')

#make_cluster_report_csv(dr+'/'+'ncbi_plus_10X_Bo.E.Heads.pass.corrected.stranded.choped.modified.cutadapt_tofu.collapsed.group.txt')
#make_collapsed_abundance_txt(dr+'/'+'ncbi_plus_10X_Bo.E.Heads.pass.corrected.stranded.choped.modified.cutadapt_tofu.collapsed.group.txt',dr+'/'+'Bo.E.Heads.pass.corrected.stranded.choped.modified.cutadapt.fasta')
#make_collapsed_abundance_txt2(dr+'/'+'all_samples.chained_count.txt', 21442936)
#make_collapsed_abundance_txt2(dr+'/'+'CAJHJT01_C1-15H_pass_FL_canu.cutadapt.unmapped.flair.NEW.firstpass.q.counts_ed', 21442936)
#man_gft(dr+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gffread_delete.gtf',dr2+'/'+'coverage.1H.bed')

#man_gft(dr+'/'+'1234')
#man_gft2('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/unmapped/novelgenes/sqanti/EGII_novelgenes.tofu.collapsed.flair.sqanti_corrected.gtf.cds.gff')
#make_fasta_single_line(dr+'/'+'protein.gbk.fasta')
#make_fasta_single_line(dr+'/'+'single_protein.gbk.fasta.sorted')
#add_uniprot_id_to_clusters(dr+'/'+'b_oleae_expression_over_1500_out_optimal_clustering.txt','/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/ncbi_genome/full_annotation/proteins/functional_annotation/Bo.all_proteins_UP000000803_blastp_sorted.txt')
#check_bad_gtf(dr+'/'+'GCF_000347755.2_Ccap_1.1_genomic.gtf')
#add_stuff(dr+'/'+'dmel.proteins.lengths',dr+'/'+'Bo.all_proteins_UP000000803_blastp_sorted.txt')

#get_genomic_ranges(dr+'/'+'ncbi_importins_blast.sort.txt', dr+'/'+'readlengths')
#get_genomic_ranges(dr+'/'+'10X_74X_importins_blast.sort.txt', dr+'/'+'readlengths')
#get_genomic_ranges(dr+'/'+'10X_92X_importins_blast.sort.txt', dr+'/'+'readlengths')
#combine_lines(dr+'/'+'mapped_reads_refs')
#bed12togtf(dr+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds.gff3.mrna.gencode.gtf.collapsed.genePred.sorted.introns.bed12')

#find_man_identifier(dr+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds.rna.fa')
#alignqc_bed(dr+'/'+'chimera.bed')
#file_manip_u_fool(dr)
#collapse_overlaps('/home/banthony/scratch/analysis/combined/10x_All/polished_50X/miRNA/bowtie/1')
#collapse_overlaps('/home/banthony/scratch/analysis/combined/10x_All/polished_50X/miRNA/bowtie/MU_boleae_V2_insects.mature_miRNA_with_alignments.sorted.bed')
#cant_use_this_again('/home/banthony/scratch/analysis/combined/10x_All/polished_50X/orthologs/orthofinder/OrthoFinder/Results_Jul09_1/Orthogroups/Orthogroups.GeneCount.tsv')
#inutile_xyz('Boleae_MU_v2_TEannot.final.gff3')
#inutile_xyz('zGalaxy304-Run1-TEannot.gff3_edited')
#nevertouseagain('cleaned.fasta')
#find_low_trembl(dr+'/'+'Bo.all_proteins_UP000000803_blastp.txt')
#get_genomic_ranges2('/home/banthony/scratch/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/gfold/maternal_degraded_completely', 
					#'/home/banthony/scratch/analysis/nanopore/olivefly/ncbi_genome/full_annotation/ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gff3',
					#'/home/banthony/scratch/analysis/nanopore/olivefly/ncbi_genome/full_annotation/ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs_10X_only_gene_refs.fasta')
#pick_random('/home/banthony/scratch/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/gfold/maternal_degraded_zygotic')
###man_gft3(dr+'/'+'EGII_novel_genes.tofu.collapsed.gff')
#see_how_bad('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/newgtf/3')
#reverse_some_seqs('/home/banthony/scratch/analysis/covid19/entebbe_primers/RPM_B_selected_primers.fasta')
#reverse_some_seqs('/home/banthony/scratch/analysis/covid19/entebbe_primers/artic_V3_reverse.fasta')
#reverse_some_seqs('/home/banthony/scratch/analysis/covid19/entebbe_primers/FIOCRUZ_reverse.fasta')
#reverse_some_seqs('/home/banthony/scratch/analysis/covid19/entebbe_primers/netherlands_AI_reverse.fasta')
#reverse_some_seqs('/home/banthony/scratch/analysis/covid19/entebbe_primers/SARS-Cov2-Midnight-1200_reverse.fasta')
#reverse_some_seqs('/home/banthony/scratch/analysis/covid19/entebbe_primers/artic_v4_reverse.fasta')
#manipulate_sam_lines(dr+'/'+'mneonegreen_wt_class_seqs.reversed_minimap2.sam',dr+'/'+'mneonegreen_wt.fasta')
#get_longest_orfx(dr+'/'+'orf_output_sl.fasta')
#dist_reads4('/home/banthony/scratch/reads/nanopore/capitata/C011_9_15/C1H/pass/cutadapt/C1_pass_canu.FL_corNuncor_reads.cutadapt.fasta')
file1='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/new_flair/tofu/sqanti/C1-15H_pass_FL_canu.flair.sqanti_corrected.gtf.cds.gff'
file2='/home/banthony/scratch/analysis/nanopore/capitata/ncbi_genome/annotation/GCF_000347755.3_Ccap_2.1_genomic.gtf'
file2='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_FINAL_flair.gtf'
file3='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/new_flair/tofu/sqanti/C1-15H_pass_FL_canu.flair.sqanti_classification.txt'
file4='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/new_flair/tofu/C1-5H_pass_filtered_lite.flair.isoforms.tofu.collapsed.group.txt'
#file3='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/tofu/reflair/sqanti/tofu/sqanti/transcript_ids_for_flair_isoforms.tofu.collapsed.rep.sqanti_classification.txt'
#file3='/home/banthony/P37.classn'
#sqanti_gff_edit(file1,file2,file3,file4)
file1='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/C0_C15H_bulk.stringtie2.sqanti_corrected.gtf.cds.gff'
file2='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/C0_C15H_bulk.stringtie_ed.gtf'
file3='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/C0_C15H_bulk.stringtie2.sqanti_classification.txt'
#sqanti_gff_edit2(file1,file2,file3)
file1='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/tofu/reflair/sqanti/side_job/missed_transcripts_with_FLAIR_associated_genes.tomoveforward'
file2='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/tofu/reflair/sqanti/side_job/transcripts_with_FLAIR_associated_genes'
classn='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/tofu/reflair/sqanti/side_job/missed_transcripts_with_FLAIR_associated_genes.sqanti_classification.txt'
ncbi_gtf='/home/banthony/scratch/analysis/nanopore/capitata/ncbi_genome/annotation/ERCC92_GCF_000347755.3_Ccap_2.1_genomic.gtf'
sqanti_gff='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/tofu/reflair/sqanti/side_job/missed_transcripts_with_FLAIR_associated_genes.sqanti_corrected.gtf.cds.gff'
#man_gft_x4(file1, file2, classn, ncbi_gtf, sqanti_gff)
file1='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/new_flair/tofu/sqanti/C1-15H_pass_FL_canu.flair.sqanti_corrected.gtf.cds_updated.gtf'
file2='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/new_flair/tofu/sqanti/C1-15H_pass_FL_canu.flair.sqanti_corrected.fasta'
#man_gft_x5(file1,file2)
#delete_now('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/new_flair/tofu/sqanti/z1')
#delete_zxxz('/home/banthony/scratch/analysis/delete/diff','/home/banthony/scratch/analysis/delete/genes')
file1='/home/banthony/scratch/analysis/covid19/entebbe_primers/blast/artic_V3_forward_blast'
file2='/home/banthony/scratch/analysis/covid19/entebbe_primers/blast/artic_V3_reverse_blast'
file1='/home/banthony/scratch/analysis/covid19/entebbe_primers/blast/artic_V4_forward_blast'
file2='/home/banthony/scratch/analysis/covid19/entebbe_primers/blast/artic_V4_reverse_blast'
file3='/home/banthony/scratch/analysis/covid19/entebbe_primers/variants/B.1.1.7'
#find_overlapping_ranges(file1,file2,file3)
file3='/home/banthony/scratch/analysis/covid19/entebbe_primers/variants/B.1.1.529'
#find_overlapping_ranges(file1,file2,file3)
file3='/home/banthony/scratch/analysis/covid19/entebbe_primers/variants/B.1.351'
#find_overlapping_ranges(file1,file2,file3)
file3='/home/banthony/scratch/analysis/covid19/entebbe_primers/variants/B.1.617'
#find_overlapping_ranges(file1,file2,file3)
file3='/home/banthony/scratch/analysis/covid19/entebbe_primers/variants/B.1.617.2'
#find_overlapping_ranges(file1,file2,file3)
file3='/home/banthony/scratch/analysis/covid19/entebbe_primers/variants/P.1'
#find_overlapping_ranges(file1,file2,file3)
file3='/home/banthony/scratch/analysis/covid19/entebbe_primers/variants/A.2.5.2'
#find_overlapping_ranges(file1,file2,file3)
#walk_dir1_gunzip('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P3/20220516_2223_2-A9-D9_PAI40249_47ae1593/fastq_pass/')
#walk_dir1_gunzip('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6C1/20220517_2009_2-A11-D11_PAH99495_ff7a3112/fastq_pass/')
#walk_dir1_gunzip('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_9_12/Pool2/fastq_pass/')
#walk_dir2_porechop('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_9_8_P12/fastq_pass/Panh/7/')
#walk_dir2_porechop('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P1/20220516_2222_2-A1-D1_PAI53876_0827877c/fastq_pass/')
#walk_dir2_porechop('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P2/20220516_2222_2-A7-D7_PAI40217_0be44d3b/fastq_pass/')
#walk_dir2_porechop('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P3/20220516_2223_2-A9-D9_PAI40249_47ae1593/fastq_pass/')
#walk_dir3('/home/banthony/scratch/analysis/covid19/B004/', '/home/banthony/scratch/analysis/covid19/B004/barcode_assignment')
#walk_dir3_count_reads('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/fastq_fail', '/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/barcode_assignment')
#walk_dir3_count_reads('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/fastq_pass', '/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/barcode_assignment')
#return_total_reads('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/fastq_pass/', '/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/barcode_assignment')
#walk_dir3_flagstat('/home/banthony/projects/rrg-ioannisr/banthony/analysis/covid19/B004_twist1/','/home/banthony/scratch/reads/nanopore/covid19/B004/B004_twist1/barcode_assignment')
#walk_dir3_flagstat2('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/','/home/banthony/scratch/reads/nanopore/capitata/barcode_assignment_ALLplates')
#walk_dir3_flagstat3('/home/banthony/scratch/reads/nanopore/capitata/','/home/banthony/scratch/reads/nanopore/capitata/barcode_assignment_ALLplates')
#walk_dir3_flagstat4('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/','/home/banthony/scratch/reads/nanopore/capitata/barcode_assignment_ALLplates')
#manip_samfile2('/home/banthony/scratch/analysis/covid19/B004/minimap2/MT007544_allreads_minimap2.sam','/home/banthony/scratch/analysis/covid19/B004/barcode_assignment')
#move_files('/home/banthony/scratch/reads/nanopore/covid/B004/','/home/banthony/scratch/analysis/covid19/B004/B004_twist1/barcode_assignment')
#move_files_unzip('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P3/20220516_2223_2-A9-D9_PAI40249_47ae1593/fastq_pass/', '/home/banthony/scratch/reads/nanopore/capitata/barcode_assignment_plate3')
#make_adapter_file('/home/banthony/software/Porechop/porechop/SARS-CoV-2_primers')
#cigar_to_position('/home/banthony/scratch/analysis/covid19/entebbe_primers/omicron/test.sam')
#manipulate_summary_seq_file2('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P3/20220516_2223_2-A9-D9_PAI40249_47ae1593/sequencing_summary_PAI40249_4cfd44dc.txt','/home/banthony/scratch/reads/nanopore/capitata/barcode_assignment_ALLplates')
#manipulate_summary_seq_file2('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P1/20220516_2222_2-A1-D1_PAI53876_0827877c/sequencing_summary_PAI53876_b368a32b.txt','/home/banthony/scratch/reads/nanopore/capitata/barcode_assignment_ALLplates')
#manipulate_summary_seq_file2('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P2/20220516_2222_2-A7-D7_PAI40217_0be44d3b/sequencing_summary_PAI40217_aba8f5d7.txt','/home/banthony/scratch/reads/nanopore/capitata/barcode_assignment_ALLplates')
#entebbe_primers_to_bed('/home/banthony/software/fieldbioinformatics/test-data/primer-schemes/nCoV-2019/Entebbe/amplicon_details')
#run_artic('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/fastq_pass/artic_v3','/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/barcode_assignment')
#run_artic('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/fastq_pass/entebbe/','/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/barcode_assignment')
#run_artic('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/fastq_pass/snap/','/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/barcode_assignment')
#collect_metrics('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/fastq_pass/','/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/barcode_assignment')
#walk_combine_fastq('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/fastq_pass/snap')
#walk_summarise_covseq('/home/banthony/scratch/analysis/covid19/B004_9_8_P12/Panh/')
#walk_getLineage_summaries('/home/banthony/scratch/analysis/covid19/B004_9_4/summarise_omi/')
#walk_getLineage_summaries('/home/banthony/scratch/analysis/covid19/B004_13_1/summarise_qiaseq/')
#walk_summarise_snap('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/fastq_pass/snap/')
#copy_bams_rename2('/home/banthony/scratch/analysis/covid19/B004_13_1/')
#copy_bams_rename3('/home/banthony/scratch/analysis/covid19/B004_13_1/summarise_qiaseq/','/home/banthony/scratch/analysis/covid19/B004_9_4/summarise_omi/')
#correct_qiaseq('/home/banthony/software/fieldbioinformatics/qiaseq/primer-schemes/nCoV-2019/V3/QIAseqDIRECTSARSCoV2primersfinal.bed')
#quick_merge('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_3/fastq_pass2/')
#R2C2('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_3/fastq_pass/artic_v4C/')
#split_R2C2('/home/banthony/scratch/analysis/covid19/B004_9_4/preprocess/lambda_phage_splint/R2C2_raw_reads.fastq')
#walk_split_R2C2('/home/banthony/scratch/analysis/covid19/B004_9_8_P12/Panh/preprocess4/lambda_phage_splint/')
#return_read_lengths_details('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_13_1/fastq_pass/qiaseq/','/home/banthony/scratch/reads/nanopore/covid19/B004/B004_13_1/barcode_assignment2_ed')
#return_read_lengths_details2('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P1/20220516_2222_2-A1-D1_PAI53876_0827877c/fastq_pass/','/home/banthony/scratch/reads/nanopore/capitata/barcode_assignment_ALLplates')
#artic_to_custom('/home/banthony/scratch/SARS-CoV-2.primer.bed')
#lampore_analysis('/home/banthony/scratch/reads/nanopore/covid19/EPO/')
#edit_gtf('/home/banthony/projects/rrg-ioannisr/banthony/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_FINAL_flair.gtf')
#combine_liqa_results('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo')
#corinne_file_man('/home/banthony/scratch/analysis/covid19/EPO/test/corinne_test2')
#corinne_file_man2('/home/banthony/scratch/analysis/covid19/EPO/test/corinne_test2')
#corinne_file_man3('/home/banthony/scratch/analysis/covid19/EPO/test/corinne_final_file')

#return_genome_coverage2('/home/banthony/scratch/analysis/covid19/B004_11_2/plot_amplicons/1')
#return_genome_mean_cov('/home/banthony/scratch/analysis/covid19/B004_11_2/plot_amplicons/')
#return_genome_mean_cov('/home/banthony/scratch/analysis/covid19/B004_13_1/plot_amplicons/qiaseq')
#select_aligned_reads('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_twist1/B004_twist1x/')
#select_aligned_reads3('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_3/fastq_pass/')
#select_aligned_reads('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_twist1/B004_twist1x/')
#select_aligned_reads2('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/fastq_pass/')
#return_read_names('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/sra')
#return_read_names('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_twist1/sra')
#return_read_names('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_3/sra')
#combine_fastq_files('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_3/sra','/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_3/sra/Doc3')
#combine_fastq_files('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_12/sra','/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_3/sra/Doc4')
#combine_fastq_files('/home/banthony/scratch/reads/nanopore/covid19/B004/B004_twist1/sra','/home/banthony/scratch/reads/nanopore/covid19/B004/B004_11_3/sra/Doc2')
#select_unligned_reads4('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P3/20220516_2223_2-A9-D9_PAI40249_47ae1593/fastq_pass/')
#select_unligned_reads4('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P2/20220516_2222_2-A7-D7_PAI40217_0be44d3b/fastq_pass/')
#select_unligned_reads4('/home/banthony/scratch/reads/nanopore/capitata/Bulk/')
#collect_align_fastq('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P1/20220516_2222_2-A1-D1_PAI53876_0827877c/fastq_pass/','C0_C5H')
#filter_fastq_rdlen('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P1/20220516_2222_2-A1-D1_PAI53876_0827877c/fastq_pass/combined.all.porechop.pychop.cutadapt.fq')
#man_gft_x6('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/medfly_Novel_transcripts.stringtie_ed.gtf')
#man_gft_x6('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/C0_C15H_bulk.stringtie_ed.gtf')
#man_gft_x7('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_FINAL_flair.gtf')
#man_gft_x8('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_FINAL_flair.gtf_ed')
dirc='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/'
dirc2='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/'
dirc='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/sqanti_rep/'
dirc2='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/'
#man_gft_x9(dirc+'/'+'C0_C15H_bulk.stringtie.sqanti_corrected.gtf.cds.gff',dirc2+'/'+'medfly_Novel_transcripts.stringtie_ed.gtf',dirc+'/'+'Sqanti_NOVEL_transcripts.txt', dirc2+'/'+'medfly_Novel_transcripts.stringtie_ed_Novel.gtf')
#man_gft_x9(dirc+'/'+'C0_C15H_bulk.stringtie2.sqanti_corrected.gtf.cds.gff',dirc2+'/'+'C0_C15H_bulk.stringtie_ed.gtf',dirc+'/'+'Sqanti_NOVEL_transcripts.txt', dirc2+'/'+'C0_C15H_bulk.stringtie_ed_Novel.gtf')
#run_stringtie('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P1/20220516_2222_2-A1-D1_PAI53876_0827877c/fastq_pass/')
#run_stringtie('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P2/20220516_2222_2-A7-D7_PAI40217_0be44d3b/fastq_pass/')
#run_stringtie('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P3/20220516_2223_2-A9-D9_PAI40249_47ae1593/fastq_pass/')
#man_gft_x10('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_flair_tofu_stringtie.gtf')
#man_gft_x10('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_flair_tofu_stringtie.sorted.gtf')
gtf1='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_flair_tofu_stringtie.sorted.gtf'
gtf2='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/sqanti_rep/QC/filter/NOVEL_genes.stringtie.sqanti.filter.filtered.gtf'
file1='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/sqanti_rep/QC/filter/NOVEL_genes.stringtie.sqanti.filter_inclusion-list.txt_ed'
#man_gft_x11(gtf1, gtf2, file1)
#man_gft_x12(gtf1)
#man_gft_x12(gtf1.replace('.gtf','_ed.gtf'))
#man_gft_x13('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_flair_tofu_stringtie.sorted_ed.gtf')
#combine_stringtie_quantification('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie_expn/all/')
#man_gft_x14('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_flair_tofu_stringtie.sorted_ed.gtf')
rsem='/home/banthony/scratch/analysis/illumina/capitata/sra/rsem/GCF_000347755.3_Ccap_2.1_ALL_rsem.genes.results'
stringtie='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie_expn/C0_15H.combined.stringtie.tpm'
#combine_rsem_to_stringtie(rsem,stringtie)
#multiply_reads('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/NOVEL_transcripts.fa')
#flair_read_matrix('/home/banthony/scratch/reads/nanopore/capitata/B002_05_6_P3/20220516_2223_2-A9-D9_PAI40249_47ae1593/fastq_pass/','/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie_expn/sampleid_sex','batch3')
file1='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie/mapped/sqanti/final_files/sqanti_rep/QC/filter/NOVEL_genes.stringtie.sqanti.filter_RulesFilter_result_classification.txt'
file2='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/flair/counts_matrixABC_bulk_ed.tpm.tsv'
#add_average_tpm(file1,file2,0.7)
faa='/home/banthony/Documents/analysis/capitata/rna_seq/single_embryo2/Cap_all/dmel_homologs/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_flair_tofu_stringtie.sorted_ed.NovelGenes.faa'
bed='/home/banthony/Documents/analysis/capitata/rna_seq/single_embryo2/Cap_all/dmel_homologs/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_flair_tofu_stringtie.sorted_ed.NovelGenes_longest.bed'
#get_longest_orf(faa, bed)
dr = '/home/banthony/R/tests/test49_rep2/absolute_quantification/'
#reorder_columns(dr+'C0_15H.combined.stringtie.tpm.rsem.filtered.tpm', dr+'correct_order_samplenames')
#reorder_columns(dr+'counts_matrixABC_bulk_ed.flair.filtered.tsv', dr+'correct_order_samplenames')
#reorder_columns(dr+'counts_matrixABC_bulk_ed.flair.tpm.filtered.tsv', dr+'correct_order_samplenames')
dr = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie_expn/'
#change_stringtie_colnames(dr+'C0_15H.combined.stringtie.tpm.rsem.filtered.tpm_ed',dr+'correct_order_samplenames')
calc_gene_expn('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/stringtie_expn/C0-15H_combined.stringtie.tpm.filtered.TPE.tsv')