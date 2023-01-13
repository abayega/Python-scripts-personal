#!/usr/bin/env python3

import shutil, re, os, sys, shelve #, patternCount5
import subprocess
import progressbar
import check_fasta

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
		if not "NC_000857" in name:
			fastaNameSeq.setdefault( name.split()[0].strip()+'\n', sequence.strip()+'\n')
			#fastaNameSeq.setdefault(name.split()[0].split('_')[0].strip('>').strip(), sequence.strip()+'\n')
		#fastaNameSeq.setdefault(name.split()[0].strip()+'\n', sequence.strip()+'\n') #fastaNameSeq.setdefault(name.split()[0].strip()+'\n', sequence.strip()+'\n')
		#name_ed = name.strip('>').split()[0]
		#fastaNameSeq.setdefault(name_ed.strip()+'\n',sequence.strip()+'\n') #fastaNameSeq.setdefault(name_ed,name.strip()+';'+sequence.strip()+'\n')		
		'''
		aname = name.split()[0].strip()+'\n' #name.split('_')[-1] #"_".join(name.split())
		if aname not in fastaNameSeq.keys():
				fastaNameSeq.setdefault(aname,[sequence.strip()+'\n'])
		elif aname in fastaNameSeq.keys():
				fastaNameSeq[aname].append(sequence.strip()+'\n')
		'''
		'''
		try:
			nompat = re.compile(r'GB_protein:(\w+.\d+)')
			nompat = re.compile(r'FlyBase_Annotation_IDs:(\w+-\w+)')
			aname = nompat.search(name).group(1)+'\n'
			fastaNameSeq.setdefault('>'+aname, sequence.strip()+'\n')
		except:
			print(name)
		'''
		#fastaNameSeq.setdefault(aname+'\n', sequence.strip()+'\n')
		#fastaNameSeq.setdefault(aname+'\n', "_".join([name.split()[0],aname+'\n'])+sequence.strip()+'\n')		
		x+=2
		bar.update(x)
	bar.finish()
	
	#print(fastqNameSeq)
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
		a1=a.split()[0]+'\n'
		b=file1.readline()
		c=file1.readline()
		d=file1.readline()
		
		#fastqNameSeq.setdefault(a, b+c+d.strip()+'\n')
		fastqNameSeq.setdefault(a1, a+b+c+d.strip()+'\n')
		x+=4
		bar.update(x)
	bar.finish()
	
	#print(fastqNameSeq)
	return(fastqNameSeq)

def look_up_given_reads_fasta(file,fasta): #Use this for speed
	#This script takes a fasta file and a file with the names of reads u want to take from the fasta. It will create a dictionary for the fasta seqs and search read names against the fasta name. 
	#If a name in the file is found in the fasta as well, the corresponding sequence and the name are returned. Pay attention to the \n and > and fasta format, it has to be single line for seqs
	
	print("\nCreating dictionary of the fasta...")
	#lopenfasta = fasta_to_dict(fasta)
	lopenfasta = fastq_to_dict(fasta)
	openfile=open(file,'r')
	lopenfile = openfile.readlines()
	l = len(lopenfile)
	path = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/new_flair' #os.path.dirname(file); name = os.path.basename(file)
	outa=open(path+'/'+'C1-15H_pass_FL_canu.flair.sqanti_classification.filtered_lite.multi_exon.fasta','w')
	
	#outb=open(path+'/'+'uncorrected_reads'+'.fasta','w')
	#Check
	print([lopenfile[0]])
	print(list(lopenfasta.keys())[0])
	total_reads = len(lopenfasta.keys())
	#print(['>'+lopenfile[0].strip().split('\t')[0]+'\n']) #print([lopenfile[0].split(".")[0]])
	a=b=c=d=e=f=n = 0
	print("\nQuerrying fasta dictionary for reads...\n")
	yyy=isoforms=lengths = []
	ref_Pat = re.compile(r'(_\d+)') #ref_Pat = re.compile(r'_NW_(\d+.\d+:\d+)')
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for read in lopenfile:
		#gene = read.strip().split()[0]
		#read = read.split(".")[0]+'\n'
		#read = read.split()[0]+'\n'
		#read = '>'+read.strip().split()[1]+'\n'
		#read = "_".join(read.strip().split()) #read.strip().split()[0]+'_'+read.strip().split()[1] +'\n'

		readnamex = '>'+read.split()[0].strip()+'\n'
		#if 'rna' not in readnamex and 'gene' not in readnamex and readnamex in lopenfasta.keys():
			#refx = ref_Pat.search(readnamex).group(1)
			#readnamey = readnamex.replace(refx,'')
			#a += 1; outa.write('>'+read.split()[0]+'\n' + lopenfasta[readnamex]) #outa.write(readnamex + lopenfasta[readnamex]) #outa.write('>'+readid+sequence)
			
		if readnamex in lopenfasta.keys(): #'>'+read.strip().split()[0]+'\n' in lopenfasta.keys():
			outa.write(lopenfasta[readnamex])
			#read_ed = lopenfasta[read.strip()].split(';')[0]
			#sampleid = read_ed.split()[2].replace('sampleid=','')
			#sequence = lopenfasta[read.strip()].split(';')[1] #geneTrans = read.strip()+'_'+gene+'\n'
			#outa.write('>'+read.strip()+'_'+sampleid+'_b4corr\n'+sequence) #print(read+lopenfasta[read].replace('*',''))
			#lengths += [str(len(sequence.strip()))+'\n']
			#outa.write(str(geneTrans)+str(lopenfasta[read].replace('*','')))
			#print(str(read)+str(lopenfasta[read].replace('*','')).strip())
			#del lopenfasta['>'+read.strip().split()[0]+'\n']
			a += 1
		'''
		meta = read.strip().split('\t')
		Genebank = read.strip().split('\t')[3].split(',')
		
		for xp in Genebank:
			if '>'+xp+'\n' in lopenfasta.keys():
				read_name = '>'+'|'.join([meta[0],meta[1],str(len(lopenfasta['>'+xp+'\n'].strip())),xp,meta[2]])+'\n'
				outa.write(read_name+lopenfasta['>'+xp+'\n'].strip()+'\n')
				yyy += ['>'+xp+'\n']
				x += 1
			elif '>'+xp+'\n' not in lopenfasta.keys() and not xp == '.':
				print(xp)
		'''
		'''
		for name,seq in lopenfasta.items():
			if read in name:
				outa.write(name+seq)
				x += 1
		'''
		'''
		if read in lopenfasta.keys():
			isoforms += [len(lopenfasta[read])]
			outa.write(''.join(lopenfasta[read]))
			x += 1
		'''
		'''
		if read in lopenfasta.keys():
			del lopenfasta[read]
			x += 1
		'''
		n+=1
		bar.update(n)
	bar.finish()
	#outc=open(path+'/'+'readLengthB4Corr','w')
	#outc.write(''.join(lengths))
	#outc.close
	#lengths = []
	'''
	rr=0
	for name,read in lopenfasta.items():
		if name.startswith('>i'):
			outa.write(name+read)
		#readname = item.split(';')[0]
		#sampleid = readname.split()[2].replace('sampleid=','')
		#sequence = item.split(';')[1]
		#lengths += [str(len(sequence.strip()))+'\n']
		#outb.write(readname.split()[0]+'_uncorr\n'+sequence)
			b += 1
		else:
			c += 1
	'''	
		
		#if name == '>Contig3811_pilon_pilon\n':
		#	seq = seq[30861:len(seq)]
		#	outa.write(name+seq)
		#elif name == '>Contig2478_pilon_pilon\n':
		#	seq = seq.strip()[0:785407].strip('N')+'\n'
		#	outa.write(name+seq)		
		#else:
		#	outa.write(name+seq)	
	#'''
	#outc=open(path+'/'+'readLengthUncorr','w')
	#outc.write(''.join(lengths))
	#outc.close
	#outa.close()
	outa.close()
	print('Nice reads = ' + str(b) + '\nSingle exon reads = ' + str(a) + '\nMatched reads = ' + str(c)) # + '\nUncorrected reads = ' + str(b) + '\nCorrected reads not found= '+str(l-a))
	#print('Percentage corrected reads = ' + str(round(100*a/l,0)))
	#print('\nDiscovered ' + str(x) + ' sequences' + ' which is %s percent of all the names u gave' %(str(round(100*x/l))))
	#print(sorted(isoforms))

def look_up_given_reads_fasta4(file, file2, fasta, fasta2): #Use this for speed
	#This script takes a fasta file and a file with the names of reads u want to take from the fasta. It will create a dictionary for the fasta seqs and search read names against the fasta name. 
	#If a name in the file is found in the fasta as well, the corresponding sequence and the name are returned. Pay attention to the \n and > and fasta format, it has to be single line for seqs
	
	print("\nCreating dictionary of the fasta...")
	lopenfasta = fasta_to_dict(fasta)
	lopenfasta2 = fasta_to_dict(fasta2)
	openfile2=open(file2,'r')
	lopenfile2 = openfile2.readlines()
	lopenfile2_dic = {}
	for linex in lopenfile2[1:]:
		lopenfile2_dic.setdefault(linex.strip().split()[0],float(linex.strip().split()[5]))
	openfile=open(file,'r')
	lopenfile = openfile.readlines()
	l = len(lopenfile)
	path = os.path.dirname(fasta)
	name = os.path.basename(fasta)
	outa=open(path+'/'+'flair_FINAL_sqanti_plus_manually_filtered_isoforms'+'.fasta','w') #; outb=open(path+'/'+'flair_exactMatches'+'.fasta','w'); outc=open(path+'/'+'flair_new_isoforms'+'.fasta','w') #outa=open(path+'/'+name+'.fasta','w')
	
	#Check
	print([lopenfile[0]])
	print([list(lopenfasta.keys())[0]])
	total_reads = len(lopenfasta.keys())
	a=b=c=d=e=f=g=h=n = 0
	print("\nQuerrying fasta dictionary for reads...\n")
	isoforms=isoforms2=lengths = []
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for read in lopenfile[1:]:
		fields = read.split()
		readname = fields[0].replace('%7C','|')
		if readname.startswith('rna'): #len(readname.split('_')) == 2 or 
			a += 1
			continue
		elif 'gene' in readname: #'rna' in readname or 
			b += 1
			isoforms += [readname]
		elif float(fields[25]) >= 1:
			c += 1; isoforms2 += [readname]
		else:
			d += 1
		'''
		if read.strip().split()[0].isdecimal() and read.strip().split()[0] in lopenfile2_dic.keys():
			#print('yes ' + lopenfile2_dic[read.strip().split()[0]])
			if lopenfile2_dic[read.strip().split()[0]]+'\n' in lopenfasta.keys():
				#print('yes2')
				del lopenfasta[ lopenfile2_dic[read.strip().split()[0]]+'\n' ]
				a += 1
		'''
		n+=1
		bar.update(n)
	bar.finish()
	#print(isoforms)
	#print('>'+isoforms[0]+'\n' in lopenfasta.keys())
	for readnamex in isoforms2: #readid,sequence in lopenfasta.items():
		readnamex = '>'+readnamex+'\n'
		if not "dup" in readnamex and readnamex not in lopenfasta2.keys():
			e += 1; outa.write(readnamex + lopenfasta[readnamex][0]) #outa.write('>'+readid+sequence)
		elif "dup" in readnamex:
			dup = readnamex.strip().split('_')[-1]
			if readnamex.replace('_'+dup,'') not in lopenfasta2.keys():
				e += 1; outa.write(readnamex+lopenfasta[readnamex.replace('_'+dup,'')][int(dup.replace('dup',''))-1])
	for readid,sequence in lopenfasta2.items():
		if readid.startswith('>rna'): #len(readid.split('_')) == 2 or 
			f += 1 #; outb.write(readid+sequence[0]) #continue
		elif not readid.startswith('>rna'): #len(readid.split('_')) != 2 and #
			g += 1 #; outa.write(readid+sequence[0])
		else:
			h += 1
	outa.close() #; outb.close(); outc.close()
	print('Total reads = ' + str(total_reads) + '\nExact matches = ' + str(a) + '\nIsoforms = ' + str(b) + '\nTPM = ' +str(c) + '\nOthers not taken = ' +str(d) )
	print('Percentage bad reads = ' + str(round(100*d/total_reads,4))); print([a,b,c,d,e,f,g,h])
	
def add_transID(gff3, fasta):
	dirc = os.path.dirname(fasta)
	nfile1 = open(dirc+'/'+'GCF_000347755.3_Ccap_2.1_rna_edited.fna', 'w')
	nfile2 = open(dirc+'/'+'zerrors', 'w')
	c=d=e=f=g = 0
	items = {}
	length = int(os.popen('wc -l %s' %(gff3)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	file1=open(gff3,'r')
	lfile1=file1.readlines()
	field_3 = ['mRNA','lnc_RNA','rRNA','snoRNA','snRNA','transcript']
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line_x in lfile1:
		ftrs = line_x.split('\t')
		
		if len(ftrs) > 3 and ftrs[2] in field_3:
			rnaIDPat = re.compile(r'ID=(\w+);')
			gene_IDPat = re.compile(r'Parent=(\w+);')
			Genbank_IDPat = re.compile(r'Genbank:(\w+.\w?);')
			trans_IDPat = re.compile(r'transcript_id=(\w+.\w?)')
			
			rnaID = gene_ID = Genbank_ID = trans_ID = "."
			try:
				rnaID = rnaIDPat.search(line_x).group(1)
				gene_ID = gene_IDPat.search(line_x).group(1)
				Genbank_ID = Genbank_IDPat.search(line_x).group(1)
				trans_ID = trans_IDPat.search(line_x).group(1)
			except:
				nfile2.write(line_x)
			if 	Genbank_ID not in items.keys():
				if Genbank_ID == trans_ID:
					items.setdefault(Genbank_ID, rnaID+'_'+gene_ID+' = ')
				else:
					items.setdefault(Genbank_ID, rnaID+'_'+gene_ID+' '+trans_ID+' ')
		d+=1
		bar.update(d)
	bar.finish()
	
	lopenfasta = fasta_to_dict(fasta)
	print([len(lopenfasta.keys()),len(items.keys())])
	for Genbank,itms in items.items():
		if Genbank in lopenfasta.keys():
			seq_title = '>' + items[Genbank] + lopenfasta[Genbank]
			nfile1.write(seq_title)
			del lopenfasta[Genbank]
		else:
			nfile2.write(Genbank + ' Not found in sequences\n')
	nfile2.write('These are sequences for which we have no Genbank\n')
	nfile2.write(str(lopenfasta.keys()))
	nfile1.close(); nfile2.close()
	
def add_transID2(gtf, fasta):
	gtf_ftr = {}
	gene_IDPat = re.compile(r'gene_id "([^"]+)"')  #re.compile(r'gene_id "(.*?)"')
	trans_IDPat = re.compile(r'transcript_id "([^"]+)"')
	
	for line in open(gtf, 'r'):
		gene_ID = gene_IDPat.search(line).group(1)
		trans_ID = trans_IDPat.search(line).group(1)
		if trans_ID not in gtf_ftr.keys():
			gtf_ftr.setdefault(trans_ID,'>'+trans_ID+'_'+gene_ID+'\n')
	lopenfasta = fasta_to_dict(fasta)
	dirc = os.path.dirname(fasta)
	nfile1 = open(dirc+'/'+'GCF_000347755.3_Ccap_2.1_genomic.transcripts_flair_FINAL_edited.fasta', 'w')
	nfile2 = open(dirc+'/'+'zerrors', 'w')
	for name,seq in lopenfasta.items():
		if name in gtf_ftr.keys():
			nfile1.write(gtf_ftr[name]+seq)
		elif 'ERCC' in name:
			for x,y in gtf_ftr.items():
				if y.split('_')[1].strip() == name:
					nfile1.write(y+seq)
		else:
			nfile2.write(name+'\n')
	nfile1.close(); nfile2.close()
	print([len(gtf_ftr.keys()),len(lopenfasta.keys())])
	
def look_up_given_reads_fastq(file, fastq): #Use this for speed
	#This script takes a fastQ file and a file with the names of reads u want to take from the fasta. It will create a dictionary for the fasta seqs and search read names against the fasta name. 
	#If a name in the file is found in the fasta as well, the corresponding sequence and the name are returned. Pay attention to the \n and > and fasta format, it has to be single line for seqs
	lopenfastq = fastq_to_dict(fastq)
	openfile=open(file,'r')
	lopenfile = openfile.readlines()
	l = len(lopenfile) #len(lopenfastq.keys()) #
	path = os.path.dirname(file)
	name = os.path.basename(file)
	outq=open(path+'/'+name+'.fastq','w') #open(path+'/'+name.replace('intronic','exonic')+'.fastq','w')
	n=x = 0
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	#'''
	for read in lopenfile:
		read = '@'+read.strip()+'\n' #read.replace('_Bo_E_1','')
		#if read.replace('>','@') in lopenfastq.keys():
		#if read in lopenfastq.keys():
		#read = read.split('_')[0].replace('>','@')+'\n'
		if read in lopenfastq.keys():
			#del lopenfastq[read]
			outq.write(read+lopenfastq[read])
			#outq.write(read.replace('>','@')+lopenfastq[read.replace('>','@')])
			x += 1
		#'''
		#for xx,y in lopenfastq.items():
		#	read = xx.strip('@\n')+'_Bo_E_1\n'
		#	if not read in lopenfile:
		#		outq.write(xx+y)
		#		x += 1
		#'''
		n+=1
		bar.update(n)
	bar.finish()
	#for key,val in lopenfastq.items():
	#	outq.write(key+val)
	outq.close()
	sys.stderr.write('\nDiscovered ' + str(x) + ' sequences' + ' which is %s percent of all the names u gave \n' %(str("{0:.1f}".format(round(100*x/l),1))))
	
def negative_look_up_given_reads_fasta(file, fasta): #Use this for speed
	#This script takes a fasta file and a file with the names of reads u want to exclude from the fasta. It will create a dictionary for the fasta seqs and a list for bad reads. 
	#If a name in the file is found in the fasta as well, the corresponding sequence and the name are not returned. Pay attention to the \n and > and fasta format, it has to be single line for seqs
	a=b=c=d=0
	print("\nCreating dictionary of the fasta...")
	lopenfasta = fasta_to_dict(fasta)
	openfile=open(file,'r')
	lopenfile = lopenfile2 = openfile.readlines()
	l = len(lopenfile)
	'''
	lopenfile = []
	for line in open(file,'r'):
		lopenfile += ['.'.join(line.split('.')[:2])]
	l = len(lopenfile)
	'''
	path = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/new_flair' #os.path.dirname(file)
	name = os.path.basename(file)
	outa=open(path+'/'+'C1-15H_pass_FL_canu.flair.sqanti_corrected.multi_exon.fasta','w')
	#outb=open(path+'/'+'uncorrected_reads'+'.fasta','w')

	print([lopenfile[0]]) #print([lopenfasta.keys()[0]])
	print([list(lopenfasta.keys())[0]])
	bar = progressbar.ProgressBar(maxval=len(lopenfasta.keys()), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for name, read in lopenfasta.items():
		#if c > 100:
		#	break
		sampleid = name.split('_')[0] #.replace('sampleid=','')
		name_ed = sampleid.strip('>').split()[0] #name.split('|')[0].strip('>')+'.path1\n' #name_ed = '.'.join(name.split('|')[0].strip('>').split('.')[:2]) 
		if name_ed.strip()+'\n' not in lopenfile and name.startswith('>i'): # or name_ed not in lopenfile:
			outa.write(name+read) #outa.write(name.split('_')[0].replace('PB.','PB.FLG.')+'\n'+read)
			a += 1
		elif name_ed+'\n' in lopenfile: # or name_ed in lopenfile:
			#outa.write('>'+name_ed+'_'+sampleid+'_b4corr\n'+read)
			b += 1
			#lopenfile.remove(name_ed+'\n')
		c+=1
		bar.update(c)
	bar.finish()
	outa.close()
	#outb.close()
	print('Total reads = ' + str(c) + '\nsingle-exon = ' + str(b) + '\nNice reads = ' + str(a) + '\nCorrected reads not found= '+str(l-b))
	#print('Total reads = ' + str(b) + '\nPrinted reads = ' + str(a) + '\nBad reads = ' + str(l) + '\nUn-printed reads = '+str(b-a))
	#print(lopenfile2)
	
def longest_isoform(fasta):
	print("\nCreating dictionary of the fasta...\n")
	lopenfasta = fasta_to_dict(fasta)
	path = os.path.dirname(fasta)
	name = os.path.basename(fasta)
	suf = '.'+name.split(".")[-1]
	outa = open(path+'/'+name.replace(suf,'.longest.fasta'),'w')
	
	seqs_dic = {}
	for seqid,seq in lopenfasta.items():
		'''
		if "ORF1_" in seqid:
			seqid = seqid.strip().split()[0].split("|")[1].replace("ORF1_",'').replace(":",'\t')+'\n'
			outa.write(seqid+seq)
			if seqid.split()[1] != "0":
				print(seqid)
		'''
		#seqid = '_'.join(seqid.strip().split()[0].split(":")[0].split("_")[1:])		
		seqid = seqid.strip().split()[0].split('-')[0]
		if seqid not in seqs_dic.keys():
			seqs_dic.setdefault(seqid,[seq])
		else:
			seqs_dic[seqid].append(seq)
	for seqid,seq in seqs_dic.items():
		seq = sorted(seq, key=len)[-1]
		outa.write(seqid+'\n'+seq.replace('*',''))
	outa.close()
	
			
def look_up_given_reads_lengths(file, fastx, file_format):
	openfile=open(file,'r')
	lopenfile = openfile.readlines()
	l = len(lopenfile)
	path = os.path.dirname(file)
	name = os.path.basename(file)
	outx=open(path+'/'+name+'.lengths','w')
	
	n = 0
	x = 0	
	if file_format == 'fasta':
		lopenfasta = fasta_to_dict(fastx)
		fastxq = 0
	elif file_format == 'fastq':
		lopenfastq = fastq_to_dict(fastx)
		fastxq = 1
	else:
		print('Supply either fasta or fastq file and declare the format')
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()	
	for read in lopenfile:
		if n == 0:
			outx.write(read.strip()+'\ttrans_length\n')
		read_ed = ">"+read.split()[0]+'\n'
		if fastxq == 1 and read.replace('>','@') in lopenfastq.keys():
			outx.write(read.replace('>','@').strip()+'\t'+str(len(lopenfastq[read.replace('>','@')].split('\n')[0]))+'\n')
			x += 1
			
		elif fastxq == 0 and read_ed in lopenfasta.keys():
			#outx.write(read_ed.strip()+'\t'+str(len(lopenfasta[read_ed]))+'\n')
			outx.write(read.strip()+'\t'+str(len(lopenfasta[read_ed]))+'\n')
			x += 1
		n+=1
		bar.update(n)		
	bar.finish()
	outx.close()
	print('Discovered ' + str(x) + ' sequences' + ' which is %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*x/l),1))))
	print('Done, bye')
	
def return_read_names_lengths(fastx, file_format):
	path = os.path.dirname(fastx)
	name = os.path.basename(fastx).replace('.'+file_format, '')
	outx=open(path+'/'+name+'_readlengths.txt','w')
	
	n = 0
	xi = 0
	
	if file_format == 'fasta':
		print(file_format)
		lopenfasta = fasta_to_dict(fastx)
		l = len(lopenfasta.keys())
		fastxq = 0
		bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		bar.start()
		for x,y in lopenfasta.items():
			#outx.write(x.strip().strip(">")+'\t'+str(len(y.strip()))+'\n')
			outx.write(x.strip().strip(">").split(' ')[0]+'\t'+str(len(y.strip()))+'\n') #code added to remove > and only take first field
			#outx.write(x+':'+str(len(y.strip()))+'\n'+y)
			#outx.write(x.strip().strip(">").split(' ')[0]+'\t'+'0\t'+str(len(y.strip()))+'\n') #code added to remove > and only take first field and return bed format ranges
			#outx.write(len(y.strip())+'\n')
			xi += 1
			n+=1
			bar.update(n)
	elif file_format == 'fastq':
		print(file_format+'\t'+os.path.basename(fastx))
		lopenfastq = fastq_to_dict(fastx)
		l = len(lopenfastq.keys())
		fastxq = 1
		bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		bar.start()
		for x,y in lopenfastq.items():
			#outx.write(x.strip()+'\t'+str(len(y.split('\n')[0].strip()))+'\n')
			outx.write(str(len(y.split('\n')[0].strip()))+'\n')
			xi += 1
			n+=1
			bar.update(n)		
	bar.finish()
	outx.close()
	print('Discovered ' + str(xi) + ' sequences' + ' which is %s percent of all sequences in the file' %(str("{0:.1f}".format(round(100*xi/l),1))))
	
def return_read_names_lengths_GC(fastx, file_format):
	path = os.path.dirname(fastx)
	name = os.path.basename(fastx).replace('.'+file_format, '')
	outx=open(path+'/'+name+'_readlengths_GC.txt','w')
	
	n = 0
	xi = 0
	
	if file_format == 'fasta':
		print(file_format)
		lopenfasta = fasta_to_dict(fastx)
		l = len(lopenfasta.keys())
		fastxq = 0
		bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		bar.start()
		for x,y in lopenfasta.items():
			#outx.write(x.strip()+'\t'+str(len(y.strip()))+'\n')
			seq = y.upper().strip()
			C = seq.count("C")
			G = seq.count("G")
			A = seq.count("A")
			T = seq.count("T")
			#GC = round(100*(G + C)/len(seq))
			GC = round(100*(G + C)/(G+C+A+T))
			#outx.write(x.strip('>\n')+'\t'+str(len(seq))+'\t'+str(GC)+'\n')
			outx.write(str(len(seq))+'\t'+str(GC)+'\n')
			xi += 1
			n+=1
			bar.update(n)
	elif file_format == 'fastq':
		print(file_format+'\t'+os.path.basename(fastx))
		lopenfastq = fastq_to_dict(fastx)
		l = len(lopenfastq.keys())
		fastxq = 1
		bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		bar.start()
		for x,y in lopenfastq.items():
			#outx.write(x.strip()+'\t'+str(len(y.split('\n')[0].strip()))+'\n')
			seq = y.upper().split('\n')[0].strip()
			C = seq.count("C")
			G = seq.count("G")
			A = seq.count("A")
			T = seq.count("T")
			#GC = round(100*(G + C)/len(seq))
			GC = round(100*(G + C)/(G+C+A+T))
			outx.write(x.strip('>\n')+'\t'+str(len(seq))+'\t'+str(GC)+'\n')
			xi += 1
			n+=1
			bar.update(n)		
	bar.finish()
	outx.close()
	print('Discovered ' + str(xi) + ' sequences' + ' which is %s percent of all sequences in the file' %(str("{0:.1f}".format(round(100*xi/l),1))))

def return_read_names_lengths_GC_ratios(fastx, file_format):
	path = os.path.dirname(fastx)
	name = os.path.basename(fastx).replace('.'+file_format, '')
	outx=open(path+'/'+name+'_readlengths_GC.txt','w')
	
	try:
		fastx = check_fasta.check_fasta_fmt(fastx)
	except:
		if not os.path.exists(fastx):
			sys.exit(-1)
	n = 0
	xi = 0
	
	if file_format == 'fasta2':
		print('\n'+file_format+'\n')
		lopenfasta = fasta_to_dict(fastx)
		l = len(lopenfasta.keys())
		fastxq = 0
		bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		bar.start()
		for x,y in lopenfasta.items():
			#outx.write(x.strip()+'\t'+str(len(y.strip()))+'\n')
			#x = x.split("(LOC")[1].split("),")[0]
			seq = y.strip()
			C = seq.count("C")
			G = seq.count("G")
			GC = round(100*(G + C)/len(seq))
			outx.write(x.strip('>\n')+'\t'+str(len(seq))+'\t'+str(GC)+'\n')
			#outx.write(str(len(seq))+'\t'+str(GC)+'\n')
			xi += 1
			n+=1
			bar.update(n)
	outx.close()
	print('\nDiscovered ' + str(xi) + ' sequences' + ' which is %s percent of all sequences in the file' %(str("{0:.1f}".format(round(100*xi/l),1))))

def remove_fasta_seqs(fasta): #File must end in .fasta
	#This script takes a fasta file and it will create a dictionary for the fasta seqs. It searches the sequence names and if the sequence name has/has not a particular partern  it is printed
	#chech fasta format
	a=b=c=d=e=f=g=h=j=k=l=n=m=x=u = 0
	try:
		fasta = check_fasta.check_fasta_fmt(fasta)
	except:
		if not os.path.exists(fasta):
			sys.exit(-1)
	
	print("\nCreating dictionary of the fasta...\n")
	lopenfasta = fasta_to_dict(fasta)

	path = os.path.dirname(fasta)
	fname = os.path.basename(fasta)
	if fname.endswith('fasta'):		
		outa_clean=open(path+'/'+fname.replace('.fasta','_no_ERCC.fasta'),'w')
		outa_ercc=open(path+'/'+fname.replace('.fasta','_ERCC.fasta'),'w')
	elif fname.endswith('fa'):
		outa_clean=open(path+'/'+fname.replace('.fa','_no_ERCC.fa'),'w')
		outa_ercc=open(path+'/'+fname.replace('.fa','_ERCC.fa'),'w')
	else:
		outa_clean=open(path+'/'+fname+'_no_ERCC','w')
		outa_ercc=open(path+'/'+fname+'_ERCC','w')
	l = len(lopenfasta.keys())
	
	#Check, print([lopenfasta.keys()[0]])

	for name,read in lopenfasta.items():
		if u < 1:
			print([name])
		else:
			break
		u += 1

	print("Querrying fasta dictionary for reads...\n")
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for name,read in lopenfasta.items():
		name = name.strip()+'\n'
		read = read.strip()+'\n'
		if "ERCC" not in name:
			outa_clean.write(name+read)
			a += 1
		elif "ERCC" in name:
			outa_ercc.write(name+read)
			b += 1
		c+=1
		bar.update(c)
	bar.finish()
	outa_clean.close()
	outa_ercc.close()
	print('Discovered ' + str(a) + ' non ERCC sequences and ' + str(b) + ' ERCC sequences, the total giving %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*(a+b)/l),1))))
	
def filter_for_fasta_seqs(fasta):
	try:
		fastax = check_fasta.check_fasta_fmt(fasta)
	except:
		if not os.path.exists(fastax):
			sys.exit(-1)
	print("\nCreating dictionary of the fasta...\n")
	lopenfasta = fasta_to_dict(fastax)
	l = len(lopenfasta.keys())
	path = os.path.dirname(fasta)
	name = os.path.basename(fasta)
	outa=open(path+'/'+'tune_reads.fasta','w')
	u=n=x = 0
	#print([lopenfasta.keys()][0][0])

	print("\nQuerrying fasta dictionary for reads...\n")

	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for header,seq in lopenfasta.items():
		if 'tune' in header:
			outa.write('>'+header+seq)
			x += 1
		n+=1
		bar.update(n)
	bar.finish()
	outa.close()
	print('\nDiscovered ' + str(x) + ' sequences' + ' which is %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*x/l),2))))
	
def find_genes(gene_ids, gtf): 
	#This function takes a gtf file from sqanti and a file with gene ids. It then loops over the gtf lines checking if gene ids are in the file and printing them
	openfile=open(gene_ids,'r')
	lopenfile = openfile.readlines()
	l = int(os.popen('wc -l %s' %(gtf)).read().split()[0])
	
	path = os.path.dirname(gene_ids)
	name = os.path.basename(gene_ids)
	outg=open(path+'/'+name+'.bed','w')	
	
	n = 0
	x = 0
	
	print('\nWe are looping over ' + str(l) + ' lines\n')
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(gtf):
		'''##for sqanti gtf
		gene_id = line.strip().split()[-1].strip('";')+'\n'
		gene_id_s = '.'.join(gene_id.strip().split('.')[:2])
		if gene_id in lopenfile:
			outg.write(line.strip()+' gn_id'+' "'+gene_id_s+'";'+'\n')
			x += 1
		'''
		##for sqanti alignQC bed
		gene_id = line.split()[3]+'\n'
		if gene_id in lopenfile:
			outg.write(line)
			x += 1
		##
		n+=1
		bar.update(n)
	bar.finish()
	outg.close()
	print('\nYou provided ' + str(len(lopenfile)) + ' gene ids')
	print('\nDiscovered ' + str(x) + ' exons' + ' which is %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*x/len(lopenfile)),1))))
	
def look_up_given_reads_fasta1(file, fasta): #to create orders from PASTEC output
	try:
		fasta = check_fasta.check_fasta_fmt(fasta)
	except:
		if not os.path.exists(fasta):
			sys.exit(-1)
	print("\nCreating dictionary of the fasta...\n")
	lopenfasta = fasta_to_dict(fasta)
	openfile=open(file,'r')
	lopenfile = openfile.readlines()
	l = len(lopenfile)
	print(lopenfile[0].split('\t'))
	path = os.path.dirname(file)
	name = os.path.basename(file)
	outb=open(path+'/'+name.replace('txt','badLTRseqs')+'.txt','w')
	u = 0
	for x,y in lopenfasta.items():
		if u < 1:
			print([x])
		else:
			break
		u += 1

	n = 0
	x = 0
	print("\nQuerrying fasta dictionary for reads...\n")
	yyy = []
	final = {}
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for read in lopenfile:
		read = read.split('\t')
		ID = '>'+read[0]+'\n'
		Class = read[4]
		Order = read[5]
		Comp = read[6]
		if Order == 'noCat':
			Supfam = Order+'_'+'unknown'
			header = ":".join([ID.strip(),'Class_'+Class,Order,Supfam,Comp])+'\n'
			#header = ":".join([ID.strip(),'Class_'+Class,Order,Comp])+'\n'
		elif Order == 'LTR':
			for it in read[-1].split(';'):
				if it.strip().startswith('coding'):
					try:
						Supfam = Order+'_'+it.split()[1].split(":")[3]
					except:
						Supfam = Order+'_'+'unknown'
						#Supfam = it.split()[1].split("_")[1]+'_'+it.split()[1].split("_")[2]
						outb.write(str(it)+'\n'+Supfam+'\n'+'\n')
			header = ":".join([ID.strip(),'Class_'+Class,Order,Supfam,Comp])+'\n'
		elif Order == 'TIR' or Order == 'LINE':
			listx = [x for x in read[-1].split(';') if x.strip().startswith('coding=')]
			if len(listx) == 1:
				try:
					Supfam = Order+'_'+listx[0].split()[1].split(":")[3]
				except:
					Supfam = Order+'_'+'unknown'
					outb.write(str(read)+'\n'+'\n')
			else:
				Supfam = Order+'_'+'unknown'
			header = ":".join([ID.strip(),'Class_'+Class,Order,Supfam,Comp])+'\n'
		else:
			header = ":".join([ID.strip(),'Class_'+Class,Order,Comp])+'\n'

		if ID in lopenfasta.keys() and Order not in final.keys():
			final.setdefault(Order,[header+lopenfasta[ID]])
		elif ID in lopenfasta.keys() and Order in final.keys():
			final[Order].append(header+lopenfasta[ID])
		else:
			print(ID)
		n+=1
		bar.update(n)
	bar.finish()
	r = 0
	out=path+'/TE_orders2/'
	os.makedirs(out)
	for xheader,seq in final.items():
		outa=open(out+name.replace('txt',xheader)+'.fasta','a')
		for each_item in seq:
			outa.write(each_item)
			r += 1
		outa.close()
	print('\nDiscovered ' + str(r) + ' sequences' + ' which is %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*r/l),1))))
	
def look_up_given_reads_fasta2(file, fasta): #to create potentially autonomous TE superfamilies from PASTEC output
	try:
		fasta = check_fasta.check_fasta_fmt(fasta)
	except:
		if not os.path.exists(fasta):
			sys.exit(-1)
	
	print("\nCreating dictionary of the fasta...\n")
	lopenfasta = fasta_to_dict(fasta)
	openfile=open(file,'r')
	lopenfile = openfile.readlines()
	l = len(lopenfile)
	print(lopenfile[0].split('\t'))
	path = os.path.dirname(file)
	name = os.path.basename(file)
	outb=open(path+'/'+name.replace('txt','badLTRseqs')+'.txt','w')
	u = 0
	for x,y in lopenfasta.items():
		if u < 1:
			print([x])
		else:
			break
		u += 1

	n = 0
	x = 0
	print("\nQuerrying fasta dictionary for reads...\n")
	yyy = []
	final = {}
	bar = progressbar.ProgressBar(maxval=l, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for read in lopenfile:
		read = read.split('\t')
		ID = '>'+read[0]+'\n'
		Class = read[4]
		Order = read[5]
		Comp = read[6]
		if Order == 'noCat':
			Supfam = Order+'_'+'unknown'
			header = ":".join([ID.strip(),'Class_'+Class,Order,Supfam,Comp])+'\n'
			#header = ":".join([ID.strip(),'Class_'+Class,Order,Comp])+'\n'
		elif Order == 'LTR':
			for it in read[-1].split(';'):
				if it.strip().startswith('coding'):
					try:
						Supfam = Order+'_'+it.split()[1].split(":")[3]
					except:
						Supfam = Order+'_'+'unknown'
						#Supfam = it.split()[1].split("_")[1]+'_'+it.split()[1].split("_")[2]
						outb.write(str(it)+'\n'+Supfam+'\n'+'\n')
			header = ":".join([ID.strip(),'Class_'+Class,Order,Supfam,Comp])+'\n'
		elif Order == 'TIR' or Order == 'LINE':
			listx = [x for x in read[-1].split(';') if x.strip().startswith('coding=')]
			if len(listx) == 1:
				try:
					Supfam = Order+'_'+listx[0].split()[1].split(":")[3]
				except:
					Supfam = Order+'_'+'unknown'
					outb.write(str(read)+'\n'+'\n')
			else:
				Supfam = Order+'_'+'unknown'
			header = ":".join([ID.strip(),'Class_'+Class,Order,Supfam,Comp])+'\n'
		else:
			continue
		#else:
		#	header = ":".join([ID.strip(),'Class_'+Class,Order,Comp])+'\n'
		width = 60
		if ID in lopenfasta.keys() and Supfam not in final.keys():
			#y =lopenfasta[ID]
			#seqn = fold_seq(str(y))
			#seqn = subprocess.check_output("echo %s | fold -w %s" %(lopenfasta[ID],width), shell=True)
			#final.setdefault(Supfam,[header+seqn])
			final.setdefault(Supfam,[header+lopenfasta[ID]])
			#final.setdefault(Order,[header+lopenfasta[ID]])
		elif ID in lopenfasta.keys() and Supfam in final.keys():
			#print(lopenfasta[ID])
			#seqn = fold_seq(str(lopenfasta[ID]))
			#seqn = subprocess.check_output("echo %s | fold -w %s" %(lopenfasta[ID],width), shell=True)
			#final[Supfam].append(header+seqn)
			final[Supfam].append(header+lopenfasta[ID])
			#final[Order].append(header+lopenfasta[ID])
		else:
			print(ID)
		n+=1
		bar.update(n)
	bar.finish()
	r = 0
	out=path+'/pot_auto_TE_superfam/'
	os.makedirs(out)
	for xheader,seq in final.items():
		outa=open(out+name.replace('txt',xheader)+'.fasta','a')
		for each_item in seq:
			outa.write(each_item)
			r += 1
		outa.close()
	print('\nDiscovered ' + str(r) + ' sequences' + ' which is %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*r/l),1))))
	
def look_up_given_reads_fasta3(dirc):
	path = dirc+'/'
	if os.path.exists(path+'/out'):
		shutil.rmtree(path+'/out')
		os.makedirs(path+'/out')
	else:
		os.makedirs(path+'/out')
	fg = []
	for x in os.listdir(dirc):
		if os.path.isfile(path+x) and not x.endswith('_sl.fasta'):
			fg += [x]
	for fastax in fg:
		try:
			fasta = check_fasta.check_fasta_fmt(dirc+'/'+fastax)
		except:
			if not os.path.exists(fasta):
				sys.exit(-1)

		print("\nCreating dictionary of the fasta...\n")
		lopenfasta = fasta_to_dict(fasta)
		name = fastax
		outa=open(path+'/out/'+name.replace('fasta','collapsed.fasta'), 'w')
		u=n=r = 0
		for x,y in lopenfasta.items():
			if u < 1:
				print([x])
			else:
				break
			u += 1
		print("\nQuerrying fasta dictionary for reads...\n")

		final = {}
		bar = progressbar.ProgressBar(maxval=len(lopenfasta.keys()), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		bar.start()
		for name,seq in lopenfasta.items():
			MCL = name.split('_')[1]
			if MCL not in final.keys():
				final.setdefault(MCL,[name+';'+seq])
			elif MCL in final.keys():
				final[MCL].append(name+';'+seq)
			else:
				print(MCL)
			n+=1
			bar.update(n)
		bar.finish()
		final2 = {}
		longest = ''
		to_record = ''

		for xheader,seq in final.items():
			for each_item in seq:
				seq_name = each_item.split(';')[0]
				seq_x = each_item.split(';')[1]				
				if len(seq_x.upper().replace('N','')) > len(longest):
					to_record = each_item
					longest = seq_x.upper().replace('N','')
				r += 1
			final2.setdefault(to_record.split(';')[0],to_record.split(';')[1])
			to_record = ''
			longest = ''
		for name2,seq2 in final2.items():
			outa.write(name2+seq2.strip().strip('Nn')+'\n')
		outa.close()
		#print('\nDiscovered ' + str(r) + ' sequences' + ' which is %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*r/l),1))))
def file_manip_u_fool(file1):
	dirc = os.path.dirname(file1)
	name = os.path.basename(file1)
	nfile1 = open(dirc+'/'+name + '.new', 'w')
	'''
	for line in open(file1, 'r'):
		if line.startswith(">"):
			nfile1.write("_".join(line.strip().split())+'\n')
		elif not line.startswith(">"):
			line = line.strip().strip("Nn")+'\n'
			nfile1.write(line)
	'''
	lopenfasta = fasta_to_dict(file1)
	for seq_name,seq in lopenfasta.items():
		if seq.strip().strip("Nn") != '':
			nfile1.write(seq_name+seq)
		else:
			print(seq_name.strip()+' has no sequence and has been thrown out.')
	nfile1.close()
	
def find_matching_reads(fasta1, fasta2): #Use this for speed
	
	print("\nCreating dictionary of the fasta...\n")
	lopenfasta1 = fasta_to_dict(fasta1)
	lopenfasta2 = fasta_to_dict(fasta2)
	
	maxi = max(lopenfasta1.keys(),lopenfasta2.keys())
	l = len(maxi)
	path = os.path.dirname(fasta2)
	name = os.path.basename(fasta2)
	outa=open(path+'/'+name+'.matches','w')
	outb=open(path+'/'+name+'Errors','w')
	n = 0
	x = 0
	print("\nQuerrying fasta dictionary for reads...\n")

	#bar = progressbar.ProgressBar(maxval=maxi, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	#bar.start()
	for read,seq in lopenfasta2.items():
		try:
			name = list(lopenfasta1.keys())[list(lopenfasta1.values()).index(seq)]
			outa.write('\t'.join([read.strip(">\n"),str(len(seq)),name.strip(">")]))
			del lopenfasta1[name]
			x+=1
		except:
			outb.write(read)
			n+=1
		#if seq in lopenfasta1.values():
		#	name = list(lopenfasta1.keys())[list(lopenfasta1.values()).index(16)]
		#x+=1
		#bar.update(x)
	#bar.finish()
	outb.write(''.join(lopenfasta1.keys()))
	outa.close()
	outb.close()
	print(str(len(lopenfasta1.keys()))+'\t'+str(len(lopenfasta2.keys()))+'\t'+str(x)+'\t'+str(n))
	
def find_reads_above_phredscore(file, fileformat, score): #score is minimum. So, if u put 10, we get reads with score equal to 10 or above
	from Bio import SeqIO
	count = 0
	rec = SeqIO.parse(file,fileformat)
	rec = next(rec)
	rec.letter_annotations["phred_quality"]
	good_reads = (
		rec
		for rec in SeqIO.parse(file,fileformat)
		if min(rec.letter_annotations["phred_quality"]) >= score
	)
	count = SeqIO.write(good_reads, file.replace(fileformat, "Q%s.%s" %(str(score),fileformat)), fileformat)
	print(f"Saved {count} good reads.")
		
dr = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/tofu/reflair/sqanti/side_job'
dr2 = '/home/banthony/scratch/analysis/combined/10x_All/polished_50X/repeats/PiRATE/Final_PiRATE'
dr1 = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/unmapped/novelgenes/EGII.novelgenes.IDs'
dr1 = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/unmapped/novelgenes/EGII_novel_genes.tofu.collapsed.group.txt'
dr2 = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/unmapped/sqanti/CAJHJT01_C1-15H_pass_FL.unmapped.flair.sqanti_corrected.faa'

#look_up_given_reads_fasta(dr1,dr2)
#look_up_given_reads_fasta(dr+'/'+'corrected_read_names',dr+'/'+'100.fasta')
#look_up_given_reads_fasta(dr+'/'+'corrected_read_names',dr+'/'+'C2H_pass.fasta') #It prints out
#look_up_given_reads_fasta(dr + '/' + 'potentially_autonomous_remove', dr + '/' + 'potentially_autonomous_TEs.edited.fasta')
#look_up_given_reads_fasta(dr + '/' + 'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gffread.gtf.longest_transcript', dr + '/' + 'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.fasta2')
#look_up_given_reads_fasta('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/new_flair/tofu/mono_exon_isoforms',
#						  '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/C1-15H_pass_FL_canu.flair.sqanti_classification.filtered_lite_sl.fasta')
#look_up_given_reads_fastq(dr+'/'+'Bo-E-1H_illumina_intronic_reads.names.R1',dr+'/'+'Bo_E_1H_2-253784_S73_R1_001.cutadapt.trimmomatic.paired.fq')
#look_up_given_reads_fastq(dr+'/'+'Bo-E-1H_illumina_intronic_reads.names.R2',dr+'/'+'Bo_E_1H_2-253784_S73_R2_001.cutadapt.trimmomatic.paired.fq')
#look_up_given_reads_fastq(dr+'/'+'Bo-E-1H_illumina_exonic_reads.names.R1',dr+'/'+'Bo_E_1H_2-253784_S73_R1_001.cutadapt.trimmomatic.paired.fq')
#look_up_given_reads_fastq(dr+'/'+'Bo-E-1H_illumina_exonic_reads.names.R2',dr+'/'+'Bo_E_1H_2-253784_S73_R2_001.cutadapt.trimmomatic.paired.fq')
#dr = '/home/banthony/scratch/reads/nanopore/rna_seq/Bo_E_3H/C010_08_071117/pass'
#look_up_given_reads_fastq(dr+'/'+'Bo_E_all_pass_edited_intronicreads.Bo_3H', dr+'/'+'Bo_E_3H_C010_08_pass_edited.fastq')

#dr = '/home/banthony/scratch/reads/nanopore/rna_seq/Bo_E_4H/C010_06_191017/pass'
#look_up_given_reads_fastq(dr+'/'+'Bo_E_all_pass_edited_intronicreads.Bo_4H', dr+'/'+'Bo_E_4H_C010_06_pass_edited.fastq')

#dr = '/home/banthony/scratch/reads/nanopore/rna_seq/Bo_E_5H/combined/pass'
#look_up_given_reads_fastq(dr+'/'+'Bo_E_all_pass_edited_intronicreads.Bo_5H', dr+'/'+'Bo_E_5H_pass_edited.fastq')

#dr = '/home/banthony/scratch/reads/nanopore/rna_seq/Bo_E_6H/C010_07_4_041117/pass'
#look_up_given_reads_fastq(dr+'/'+'Bo_E_all_pass_edited_intronicreads.Bo_6H', dr+'/'+'Bo_E_6H_C010_07_pass_edited.fastq')

#dr = '/home/banthony/scratch/reads/nanopore/rna_seq/adult_female/C010_12_01_120518/pass'
#look_up_given_reads_fastq(dr+'/'+'Bo_E_all_pass_edited_intronicreads.Heads_female', dr+'/'+'Bo_Heads_female_C010_12_01_pass.fq')

dr1 = '/home/banthony/scratch/analysis/covid19/B004_9_4/process_good'
dr2 = '/home/banthony/scratch/analysis/covid19/B004_9_4/process_good'
#look_up_given_reads_fastq(dr1+'/'+'reads_to_check', dr2+'/'+'R2C2_Subreads.fastq')
dr = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo2/Cap_all/flair'
look_up_given_reads_lengths(dr + '/' + 'counts_matrixABC_bulk_ed.flair.tpm.filtered.TPE.tsv',dr + '/' + 'ERCC92_GCF_000347755.3_Ccap_2.1_genomic_flair_tofu_stringtie.fa', 'fasta')

#filter_for_fasta_seqs(dr+'/'+'boleae_nanopore_male.merged_passreads.fasta')
#return_read_names_lengths(dr+'/'+'Galaxy9-10x_MP_ONT_PacBio.0.5.pilon2_sl.noContaminants.ed2.fasta.fasta', 'fasta')
#file_manip_u_fool(dr)
#longest_isoform('/home/banthony/scratch/analysis/others/proteomes_to_compare/proteomes/Dme.fasta')
#find_matching_reads(dr1+'/'+'cleaned.fasta',dr2+'/'+'boleae_pirate_TEs_ed2.fasta')
#return_read_names_lengths_GC(dr+'/'+'fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta'
#negative_look_up_given_reads_fasta(dr+'/'+'corrected_read_names',dr+'/'+'100.fasta')
#negative_look_up_given_reads_fasta(dr+'/'+'missed_transcripts_with_FLAIR_associated_genes.tomoveforward',dr+'/'+'missed_transcripts_with_FLAIR_associated_genes.sqanti_corrected.faa')
#return_read_names_lengths_GC(dr+'/'+'fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')

#file1='/home/banthony/scratch/reads/nanopore/capitata/C011_9_14/C2H/pass/canu/C2H_pass_canu/correction/2-correction/results/0002.err'
file1='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/C1-15H_pass_FL_canu.flair.sqanti_classification.txt'
file2='/home/banthony/scratch/analysis/illumina/capitata/sra/rsem/C1H_pass_FL_canu.cutadapt.flair.isoforms.results_ed'
fasta='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/C1-15H_pass_FL_canu.flair.sqanti_corrected.fasta'
fasta2='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/C1-15H_pass_FL_canu.flair.sqanti_classification.filtered_lite.fasta'
#look_up_given_reads_fasta4(file1,file2,fasta,fasta2)
dr = '/home/banthony/scratch/reads/SRA/goke_data/'
#return_read_names_lengths_GC(dr+'/'+'barcode01/fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')
#return_read_names_lengths_GC(dr+'/'+'barcode02/fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')
#return_read_names_lengths_GC(dr+'/'+'barcode03/fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')
#return_read_names_lengths_GC(dr+'/'+'barcode04/fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')
#return_read_names_lengths_GC(dr+'/'+'barcode05/fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')
#return_read_names_lengths_GC(dr+'/'+'barcode06/fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')
#return_read_names_lengths_GC(dr+'/'+'barcode07/fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')
#return_read_names_lengths_GC(dr+'/'+'barcode08/fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')
#return_read_names_lengths_GC(dr+'/'+'barcode09/fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')
#return_read_names_lengths_GC(dr+'/'+'barcode10/fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')
#return_read_names_lengths_GC(dr+'/'+'barcode11/fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')
#return_read_names_lengths_GC(dr+'/'+'barcode12/fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed.fasta','fasta')
#return_read_names_lengths_GC(dr+'/'+'barcode45.pchop1.fastq','fastq')
#return_read_names_lengths_GC(dr+'/'+'barcode46.pchop1.fastq','fastq')
#return_read_names_lengths_GC(dr+'/'+'ERR2680377.Q10.fastq','fastq')
#return_read_names_lengths_GC(dr+'/'+'ERR2680382.Q10.fastq','fastq')
#return_read_names_lengths_GC(dr+'/'+'ERR3363660.fastq','fastq')
#return_read_names_lengths_GC(dr+'/'+'SRR8487226.Q10.fastq','fastq')
#return_read_names_lengths_GC(dr+'/'+'ERR3363658.fastq','fastq')
#return_read_names_lengths_GC(dr+'/'+'ERR6053095/ERR6053095.fastq','fastq')
							 
gff='/home/banthony/scratch/analysis/nanopore/capitata/ncbi_genome/annotation/GCF_000347755.3_Ccap_2.1_genomic.gff'
fasta='/home/banthony/scratch/analysis/nanopore/capitata/ncbi_genome/transcripts/GCF_000347755.3_Ccap_2.1_rna.fna'
#add_transID(gff,fasta)
fasta='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/newgtf/GCF_000347755.3_Ccap_2.1_genomic.transcripts_flair_FINAL.fasta'
gtf='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/newgtf/ERCC92_GCF_000347755.3_Ccap_2.1_genomic_plus_FLAIR_isoforms.gtf'
#add_transID2(gtf,fasta)
#find_reads_above_phredscore('/home/banthony/scratch/reads/SRA/goke_data/ERR2680377.fastq',"fastq",10)
#find_reads_above_phredscore('/home/banthony/scratch/reads/SRA/goke_data/ERR2680382.fastq',"fastq",10)
#find_reads_above_phredscore('/home/banthony/scratch/reads/SRA/goke_data/ERR3363660.fastq',"fastq",10)
#find_reads_above_phredscore('/home/banthony/scratch/reads/SRA/goke_data/SRR8487226.fastq',"fastq",10)
#find_reads_above_phredscore('/home/banthony/scratch/reads/SRA/goke_data/ERR6053094/ERR6053094.fastq',"fastq",10)
#find_reads_above_phredscore('/home/banthony/scratch/reads/SRA/goke_data/ERR3363658.fastq',"fastq",10)
#find_reads_above_phredscore('/home/banthony/scratch/reads/SRA/goke_data/ERR6053095/ERR6053095.fastq',"fastq",10)