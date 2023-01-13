#!/usr/bin/env python3

import re, os, sys, progressbar, check_fasta
from pathlib import Path

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
		
def find_common_entries(file1, file2):
	newfile = open(file2+'.maternal_degraded_genes', 'w')
	ofile1 = open(file1)
	lfile1 = ofile1.readlines()
	lfile11 = []
	for line in lfile1:
		lfile11 += [line.split()[0]]
	ofile2 = open(file2)
	lfile2 = ofile2.readlines()
	geneid = ''
	k=n=m= 0
	bar = progressbar.ProgressBar(maxval=len(lfile2), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in lfile2[1:]:
		#print(line.split()[8])
		#break
		if line == "" or line.startswith("#"):
			continue
		#matches = re.match('gene_id "([^"]+)"', line.split('\t')[8])
		#if matches:
		#	geneid = matches.group(1)
		#else:
		#	continue
		if line.split()[0] in lfile11: #geneid+'\n' in lfile1:
		#if not line.strip('@\n')+'_Bo_E_1\n' in lfile1: #line.split()[3]+'\n' in lfile1:
			n += 1
			#lfile2.remove(line.strip('>'))
			newfile.write(line)
			#print(len(lfile2))
		else:
			#lfile2.remove(line.strip('>'))
			k += 1
		#	print(line)
			#newfile.write(line)
		m += 1		
		bar.update(m)
	bar.finish()
	newfile.close()
	#print('There are ' + str(len(lfile1)) + ' lines in ' + file1 +
	print(str(n) + '\t' + str(k))

def find_common_entries2(file1, file2):
	newfile1 = open(file1+'_NewEd','w')
	newfile2 = open(file2+'_NewEd','w')
	gene_ids = []
	file2_dic = {}
	for line in open(file2,'r'):
		if not line.startswith('>'):
			file2_dic.setdefault(line.split()[0],line) #.split('|')[0].strip('>')+'\n')
		#file2_dic.setdefault(line.split()[0],'\t'.join([ line.strip().split()[7],line.strip().split()[9],line.strip().split()[2]])+'\n' )
		'''
		if len(line.split()) > 3 and line.split()[2] == 'exon':
			ftrs2 = line.split(";")
			#print(ftrs2)
			try:
				GeneID = [x for x in ftrs2 if x.startswith(" gene_id")]
				GeneID = GeneID[0].split('"')[1]
				gene_ids += [GeneID.strip()]
			except:
				print(ftrs2)
		'''
	print(list(file2_dic.items())[0])
	length = int(os.popen('wc -l %s' %(file1)).read().split()[0])
	x=a=b=0
	bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(file1, 'r'):
		#if x == 0:
		#	newfile.write( line.strip()+'\tnew_id\n' ) #newfile.write( line.strip()+'\texpected_count\t'+'FPKM\t'+'Fscaf\n' )
		if line.split()[0] in file2_dic.keys():
			newfile1.write( line) #.strip()+'\t'+file2_dic[line.split()[0]] ) #sys.stderr.write(line)
			newfile2.write( file2_dic[line.split()[0]])
			a += 1
		else:
			b += 1 #sys.stderr.write(line.split()[0]); sys.exit(-1)
		x += 1
		bar.update(x)
	bar.finish(); newfile1.close(); newfile2.close(); print([a,b])

def find_common_entries3(file1, file2, file3):
	dirc = os.path.dirname(file2)
	name = os.path.basename(file2)
	newfile = open(dirc+'/'+'covseq_report_edited.csv','w')
	#newfile1 = open(dirc+'/'+ 'z_missing','w')
	#newfile.write('\t'.join([ 'gene_id','length','effective_length','std_1H_in','std_1H_ex','std_2H_in','std_2H_ex','std_3H_in','std_3H_ex','std_4H_in','std_4H_ex','std_5H_in','std_5H_ex','std_6H_in','std_6H_ex','std_fem_in','std_fem_ex','std_mal_in','std_mal_ex\n' ]))
	#newfile.write('chromosome\t'+'gene_name\t'+'gene_ID\t'+'gene_length\t'+'transcript_name\t'+'transcript_length\t'+'transcript_GC\n')
	file1_dic = {}
	file2_dic = {}; file3_dic = {};irrele = []
	a=b=c=d=e=f = 0
	#'''
	for line in open(file1,'r'):
		sline = line.strip().split()
		if sline[2] not in file1_dic.keys():
			file1_dic.setdefault(sline[2],sline[3]+','+sline[4])
		elif sline[0] in file1_dic.keys():
			sys.stderr.write(line)
		a += 1
	#'''
	for liney in open(file2, 'r'):
		sline = line.strip().split(',')
		samplename = sline[0]; platform = sline[1]; GISAID_ID = sline[96]
		if samplename not in file2_dic.keys():
			file2_dic.setdefault(samplename,platform +','+ GISAID_ID
		elif samplename in file2_dic.keys():
			sys.stderr.write(line)
	basedir = '/home/banthony/scratch/analysis/covid19/B004/samples_to_analyse/'
	for liney in open(file3, 'r'):
		if liney.split(',')[0] == 'samplename':
			line1st = liney.strip()+',PANGOLIN,avgCov,grade\n'
			continue
		b += 1
		sampleID = liney.split(',')[0]
		concentration = liney.split(',')[40]
		try:
			if float(concentration) <= 31.40379: #the median
				grade = 'low'
			elif float(concentration) > 31.40379 and float(concentration) <= 84.26619:
				grade = 'mid'
			elif float(concentration) > 84.26619: #third quartile
				grade = 'high'
			c += 1
		except:
			grade = 'unknown'; d += 1
			
		if sampleID in file1_dic.keys():
			e += 1
			lineage = file1_dic[sampleID].split(',')[0]	
			to_write = liney.strip()+','+file1_dic[sampleID]+','+grade+'\n'
			if lineage not in file3_dic.keys():
				file3_dic.setdefault(lineage,{grade:[ to_write ]})
			elif lineage in file3_dic.keys():
				if grade in file3_dic[lineage].keys():
					file3_dic[lineage][grade].append(to_write)
				elif grade not in file3_dic[lineage].keys():
					file3_dic[lineage].update({grade:[ to_write ]})
	#print([file2_dic['B.1.1.7']['unknown']])
	newfile.write(line1st)
	for lin,details in file3_dic.items():		
		lineage_dir = Path(basedir,lin)
		#print([lin,lineage_dir])
		if not lineage_dir.exists():
			os.makedirs(lineage_dir)
		for gradex,lines in details.items():
			gradex_file = Path(lineage_dir, gradex+'_conc.csv')
			if not gradex_file.exists():
				gradex_filex = open(gradex_file, 'w')
				gradex_filex.write(line1st+''.join(lines))
				gradex_filex.close()
			elif gradex_file.exists() and gradex_file.is_file():
				gradex_filex = open(gradex_file, 'a')
				gradex_filex.write(line1st+''.join(lines))
				gradex_filex.close()
			newfile.write(''.join(lines))
				
	#[lin,len(file2_dic[lin['high']]),len(file2_dic[lin['mid']]),len(file2_dic[lin['low']]),len(file2_dic[lin['unknown']])]

		'''
		#low_conc = Path.join(lineage_dir, 'low_conc.csv')
		#mid_conc = Path.join(lineage_dir, 'mid_conc.csv')
		#high_conc = Path.join(lineage_dir, 'high_conc.csv')
		elif lineage_dir.exists() and lineage_dir.is_dir():
			low_conc = Path.join(lineage_dir, 'low_conc.csv')
			mid_conc = Path.join(lineage_dir, 'mid_conc.csv')
			high_conc = Path.join(lineage_dir, 'high_conc.csv')
			if not low_conc.exists():
				low = open(low_conc, 'w')
			elif low_conc.exists() and low_conc.is_file():
				low = open(low_conc, 'a')				
			if not mid_conc.exists():
				mid = open(mid_conc, 'w')
			elif mid_conc.exists() and mid_conc.is_file():
				mid = open(mid_conc, 'a')
			if not high_conc.exists():
				high = open(high_conc, 'w')
			elif high_conc.exists() and high_conc.is_file():
				high = open(high_conc, 'a')
		'''
		
	'''
	#print(list(file1_dic.items())[0])
	print('\nMoving on...\n')
	for line in open(file2,'r'):
		sline = line.strip().split() #[-1].split(';')
		if len(sline) <3 or line.startswith('#'):
			continue
		else:
			if sline[2] == 'mRNA':
				gene_IDPat = re.compile(r'Parent=(\w+)')
				gene_ID = gene_IDPat.search(line).group(1)
				GenbankPat = re.compile(r'Genbank:(\w+.\d+)')
				Genbank = GenbankPat.search(line).group(1)
				
				if gene_ID not in file2_dic.keys():
					file2_dic.setdefault(Genbank,gene_ID)
				else:
					print(gene_ID+'\t'+Genbank)
			#file2_dic.setdefault(sline[0],sline[1])
		#elif sline[0] in file2_dic.keys() and sline[2] == 'gene' or sline[0] in file2_dic.keys() and sline[2] == 'region':
			#file2_dic[sline[0]].append(sline[2]+';'+sline[-1].split(";")[0])
			#sys.stderr.write(line)
	'''
	#print(list(file2_dic.items())[0])
	'''
	for x,y in file1_dic.items():
		if x in file2_dic.keys():
			#tra = file2_dic[x].split()
			#newfile.write(y.strip()+'\t'+tra[1]+'\t'+tra[3]+'\t'+tra[4]+'\n')
			#newfile.write(y.strip()+'\t'+file2_dic[x]) #+'\n')
			newfile.write('\t'.join(y.split()[:3]+[y.split()[3]+'_Y']+y.split()[4:])+'\n')
			z += 1
		else:
			newfile.write(y+'\n')
	'''
	'''
	print('\nMoving on 2...')
	for liney in open(file2, 'r'):
		geneID = liney.split(',')[0] #genebank = liney.split()[1].strip("|").split("|")[-1]
		if geneID.split('_')[-1].startswith('dup'):
			dup = '_'+geneID.split('_')[-1]
			geneID = geneID.replace(dup,'')
		if liney.startswith('isoform'):
			newfile.write(liney.strip()+'\texpected_count\tFPKM\n')
		elif geneID in file1_dic.keys():
			newfile.write(liney.strip()+'\t'+file1_dic[geneID]+'\n') #liney.split()[0]+'\t'+file2_dic[genebank]+'\t'+genebank+'\n')
			#newfile.write(file2_dic[y.split()[0]]+'\t'+y.split()[1])
			#newfile.write(y.replace(y.split()[3],y.split()[3]+'_Y'))
			#newfile.write(file2_dic[y.split()[0]].replace('#DIV/0!','-3'))
			#newfile.write('\t'.join(y.split()[:3]+[y.split()[3]+'_Y']+y.split()[4:]))
			#for loop in file2_dic[y.split()[0]]:
			#	newfile.write(y.split()[0]+'\t'+'\t'.join(loop.split(';'))+'\t'+'\t'.join(y.split()[1].split('_'))+'\n')
			a += 1		
		else:
			newfile1.write(liney)
			b +=1
			#newfile.write(y)
	
	#print('\n'+str(100*z/len(file2_dic.keys())))
	print('\n'+str(a)+'\t'+str(b))
	'''
	newfile.close()
	print([a,b,c,d,e])
	
file1='ERCC92_GCF_000347755.3_Ccap_2.1_genomic.gtf'
file2='C1-15H_pass_FL_canu.flair.sqanti_corrected.gtf.cds.gff'
file3='flair_FINAL_sqanti_plus_manually_filtered_isoforms'
def find_common_entries4(file1, file2, file3):
	path = os.path.dirname(file1)
	newfile = open(path+'/'+ 'flair_FINAL_sqanti_plus_manually_filtered_isoforms.gtf','w')
	newfile2 = open(path+'/'+ 'zerrors','w')
	#newfile.write('chromosome\t'+'gene_name\t'+'gene_ID\t'+'gene_length\t'+'transcript_name\t'+'transcript_length\t'+'transcript_GC\n')
	file1_dic = {}
	file2_dic = {}
	z = 0
	for line in open(file1,'r'):
		gene_IDPat = re.compile(r'gene_id "([^"]+)"')  #re.compile(r'gene_id "(.*?)"')
		gene_symbPat = re.compile(r'gene_name "([^"]+)"')
		trans_symbPat = re.compile(r'transcript_id "([^"]+)"')
		trans_id = trans_symbPat.search(line).group(1)
		try:
			gene_ID = gene_IDPat.search(line).group(1)
		except:
			gene_ID = trans_id
		try:
			gene_name = gene_symbPat.search(line).group(1)
		except:
			gene_name = gene_ID
		
		if gene_ID not in file1_dic.keys():
			file1_dic.setdefault(gene_ID,gene_name)
		#elif geneID in file1_dic.keys():
		#	file1_dic[geneID].append(line.strip())
		#else:	
		#	sys.stderr.write(line)
	print('\nMoving on...\n')			
	for line in open(file2,'r'):
		trans_symbPat = re.compile(r'transcript_id "([^"]+)"')
		trans_id = trans_symbPat.search(line).group(1)		
		if trans_id not in file2_dic.keys():
			file2_dic.setdefault(trans_id,[line.strip()])
		elif trans_id in file2_dic.keys():
			file2_dic[trans_id].append(line.strip())
		else:	
			sys.stderr.write(line)
	print(list(file1_dic.items())[0])
	print(list(file2_dic.items())[0])
	print('\nNow querrying...\n')
	gene_IDPat = re.compile(r'gene_id "([^"]+)"')
	for line in open(file3,'r'):
		#id1 = line.strip().split()[0]; id2 = line.strip().split()[1]
		if line.strip() in file2_dic.keys(): # and id2 in file2_dic.keys():
			transcript_id = line.split('_')[0]
			items = file2_dic[line.strip()]
			for item_x in items:
				#newfile2.write(item_x+'\n')
				ls_item_x = item_x.split()
				gene_ID = gene_IDPat.search(item_x)
				if gene_ID:
					gene_ID = gene_ID.group(1)
				else:
					gene_ID = ''
				if gene_ID.startswith('gene') and len(gene_ID.split('_')) == 1:
					gene_name = file1_dic[gene_ID]
				elif gene_ID.startswith('gene') and len(gene_ID.split('_')) > 1:
					gene_ID = item_x.split('_')[-1].strip('";')
					if gene_ID.startswith('gene'):
						gene_name = file1_dic[gene_ID]
					else:
						gene_name = gene_ID
				elif gene_ID.startswith('novelGene_gene') :
					gene_ID = gene_ID.split('_')[1]
					gene_name = file1_dic[gene_ID]
				elif gene_ID == '':
					gene_ID = gene_name = transcript_id
				elif gene_ID.startswith('novelGene_') :
					gene_ID = gene_name = gene_ID
				else:
					gene_ID = transcript_id; gene_name = gene_ID; print(item_x)
				new_item_x = '\t'.join(ls_item_x[:7]).replace('PacBio', 'FLAIR')+'\t'+ '.\t'+ 'transcript_id "' + transcript_id + '"; ' + 'gene_id "' + gene_ID + '"; ' + 'gene_name "' + gene_name + '";'
				newfile.write(new_item_x+'\n')
			'''		
			geneid1 = item.split('_')[-1]
			gene_IDPat = re.compile(r'gene_id "([^"]+)"')  #re.compile(r'gene_id "(.*?)"')
			gene_symbPat = re.compile(r'gene_name "([^"]+)"')
			gene_id = gene_IDPat.search(line).group(1)
			gene_name = gene_symbPat.search(line).group(1)
			
			#tra = file2_dic[x].split()
			to_print = '\t'.join([file1_dic[id1],file2_dic[id2]])+'\n'
			#newfile.write(y.strip()+'\t'+tra[1]+'\t'+tra[3]+'\t'+tra[4]+'\n')
			#newfile.write(y.strip()+'\t'+file2_dic[x]) #+'\n')
			newfile.write(to_print)
			z += 1
			'''
		else:
			newfile2.write(line)
			#sys.stderr.write(line)
	newfile2.close(); newfile.close()
	print('\n'+str(100*z/len(file1_dic.keys())))
	print('\n'+str(100*z/len(file2_dic.keys())))
	print('\n'+str(z))

def find_common_entries5(file1, file2, file3, file4):
	file2_dic = {}
	full_splice_match = []
	busco_genes = []
	for line in open(file2,'r'):
		key = line.split()[0]
		file2_dic.setdefault(key,line)
	for line in open(file3,'r'):
		rna = line.strip().split()[-1]
		rna = line.strip().split()[2] #actually this is the gene
		full_splice_match += [rna]
	for line in open(file4,'r'):
		busco_genes += [line.split()[0]]
	#print(len(full_splice_match))
	dirc = os.path.dirname(file2) 
	name = os.path.basename(file1)
	nfile1 = open(dirc+'/'+name + '.maternal_degraded', 'w') #open(dirc+'/'+name + '.zygotic_early', 'w')
	nfile2 = open(dirc+'/'+name + '.maternal_degraded_gene_expn', 'w') #open(dirc+'/'+name + '.zygotic_early_gene_expn', 'w')
	length = int(os.popen('wc -l %s' %(file1)).read().split()[0])
	x=0
	bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(file1, 'r'):
		gene_id = line.strip().split()[-1]
		gene_id = line.strip().split()[3]
		'''
		if gene_id in file2_dic.keys():
			nfile1.write(line)
			nfile2.write(file2_dic[gene_id])
		x+=1
		'''
		rna_id = line.strip().split()[3]
		#rna_id = line.split()[3]
		#if gene_id in file2_dic.keys() and rna_id in full_splice_match:
		if gene_id in file2_dic.keys() and gene_id in full_splice_match and gene_id in busco_genes:
			#print(rna_id)
			#if x > 100:
			#	break
			nfile1.write(line)
			nfile2.write(file2_dic[gene_id])
		x+=1
		bar.update(x)
	bar.finish()
	nfile1.close()
	nfile2.close()
	
def find_common_entries6(file1, file2, file3, file4): #(ancestral.pub_og_id, odb9v1_OG2genes.tab, odb9v1_genes.tab, Bo.all_proteins_UP000000803_blastp_sorted.txt.plus.dmel.names)
	#I used this code to BUSCO diptera genes that correspond to our boleae dmel homologues
	file2_dic = {}
	file3_dic = {}
	file4_dic = {}
	for line in open(file2,'r'):
		key = line.split()[0]
		if not key in file2_dic.keys():
			file2_dic.setdefault(key,[line.split()[1].strip()])
		else:
			file2_dic[key].append(line.split()[1].strip())
	for line in open(file3,'r'):
		key = line.split()[0]
		if not key in file3_dic.keys():
			file3_dic.setdefault(key,[line.split()[4].strip()])
		else:
			file3_dic[key].append(line.split()[4].strip())
	for line in open(file4,'r'):
		try:
			key = line.split('\t')[15]
			if not key in file4_dic.keys():
				file4_dic.setdefault(key,[line.split()[0].split("|")[0]])
			else:
				file4_dic[key].append(line.split()[0].split("|")[0])
		except IndexError:
			break
		#	try:
		#		key = line.split('\t')[15]
		#		file4_dic.setdefault(key,line.split()[0].split("|")[0])
		#	except:
		#		print(yes)
	#print(file2_dic)
	#print(file3_dic)
	#print(file4_dic)
	dirc = os.path.dirname(file2) 
	name = os.path.basename(file1)
	nfile1 = open(dirc+'/' + 'busco_diptera_ancestral_genes', 'w')
	#nfile2 = open(dirc+'/'+ 'z.busco_diptera_ancestral_genes', 'w')
	length = int(os.popen('wc -l %s' %(file1)).read().split()[0])
	x=0
	bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(file1, 'r'):
		line = line.strip()
		if line in file2_dic.keys():
			for i in file2_dic[line]:
				value1 = i
				if value1 in file3_dic.keys():
					value2 = file3_dic[value1][0]
					if value2 in file4_dic.keys():
						for i in file4_dic[value2]:
							nfile1.write('\t'.join([i,value2,value1,line])+'\n')
		x+=1
		bar.update(x)
	bar.finish()
	nfile1.close()
	#nfile2.close()
	
def find_common_entries7(file1, file2, file3):
	dirc = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/tofu/reflair'
	newfile = open(dirc+'/master_file','w')
	newfile.write('\t'.join(['tofu_id','flair_id','scafold', 'full_length_count','norm_fl','length','effective_length','expected_count','TPM','FPKM\n']))
	#newfile.write('\t'.join(['read_name','Ref_aligned','read_length','alignment_identity','read_identity','aligned_identity','insertion','deletion',
	#						 'indels','perc_read_aligned','SAM-FLAG','aligned_string_length','chaining_score1','chaining_score2','aln_score',
	#						 'ref_edit_distance','strand_orientation','aln_type','readlength?','flag','tagname,gene,Chrpos,Chr','contig_len\n']))
	file1_dic = {}
	file2_dic = {}
	z = 0
	scafPat = re.compile(r'/\w+(N.+):'); scaf = 'gene'
	for line in open(file1,'r'):
		sline = line.strip().split()
		geneID = sline[0]
		if geneID not in file1_dic.keys() and 'transcript_id' not in line and not geneID.startswith('gene'):
			file1_dic.setdefault(geneID,'\t'.join(sline[2:7]))
		elif geneID in file1_dic.keys() and 'transcript_id' not in line and not geneID.startswith('gene'):
			sys.stderr.write(line)
	print('\nMoving on...\n')			
	for line in open(file2,'r'):
		sline = line.strip().split()
		geneID = sline[0]
		#xsome = '_'.join([sline[6],sline[7]])
		if geneID not in file2_dic.keys()  and '#' not in line and not geneID.startswith('gene'):
			file2_dic.setdefault(geneID,'\t'.join(sline[1:]))
		elif geneID in file2_dic.keys()  and '#' not in line:
			sys.stderr.write(line)
	print(list(file1_dic.items())[0])
	print(list(file2_dic.items())[0])
	print('\nNow querrying...\n')
	for line in open(file3,'r'):
		pbID = line.strip().split()[0]
		geneID = line.strip().split()[1].split(',')[0]
		if geneID.startswith('i') and 'gene' not in geneID:
			scaf = scafPat.search(geneID).group(1)
		elif geneID.startswith('i') and 'gene' in geneID:
			scaf = geneID.split('_')[-1]
		elif geneID.startswith('rna') or geneID.startswith('gene'):
			scaf = geneID.split('_')[-1]
		else:
			scaf = 'unknown'
			print(line)
				
		if geneID in file2_dic.keys() and geneID in file1_dic.keys():
			to_print = '\t'.join( [ pbID,geneID,scaf,file2_dic[geneID],file1_dic[geneID] ] )  + '\n'
			newfile.write(to_print)
			z += 1
		else:
			#newfile.write(y.strip()+'\n')
			a = 0 #; sys.stderr.write(line)
	print('\n'+str(100*z/len(file1_dic.keys())))
	print('\n'+str(100*z/len(file2_dic.keys())))
	print('\n'+str(z)); newfile.close()

def find_common_entries8(file1, file2):
	newfile = open(file2+'.bed','w')	
	file1_dic = {}
	file1_list = []
	z = 0
	for line in open(file1,'r'):
		sline = line.strip().split('\t')
		gene_IDPat = re.compile(r'gene_id "(.*?)"')
		geneID = gene_IDPat.search(line).group(1)
		if geneID not in file1_dic.keys() and len(sline) > 5 and sline[2] == 'CDS':
			file1_list += [geneID]
		'''
		geneID = sline[0]
		if geneID not in file1_dic.keys() and len(sline) > 5:
			file1_dic.setdefault(geneID+','+sline[2],sline[3])
		elif geneID in file1_dic.keys():
			sys.stderr.write(line)
		'''
	print('\nMoving on...\n')
	print('\nNow querrying...\n')
	print(file1_list[1:5])
	for line in open(file2,'r'):
		#print(line)
		'''
		id1 = line.strip().split('\t')[0]
		id2 = line.strip().split('\t')[1]
		leng = int(line.strip().split('\t')[2])
		if id1+','+id2 in file1_dic.keys():
			aln_start=int(file1_dic[id1+','+id2])
			if id2=='Contig36658_pilon_pilon':
				to_print = '\t'.join([id2,'0','2478',id1]) + '\n'
				newfile.write(to_print)
			else:
				to_print = '\t'.join([id2,str(aln_start-1000),str(aln_start+leng+1000),id1]) + '\n'
				newfile.write(to_print)
			z += 1
		else:
			#newfile.write(y.strip()+'\n')
			sys.stderr.write(line)
		'''
		if line.strip() in file1_list:
			z += 1
	#print('\n'+str(100*z/len(file1_dic.keys())))
	print('\n'+str(z))

def find_common_entries9(file1, file2, file3):
	newfile = open('boleaeGenes_to_chromosomes','w')
	newfile2 = open('investigate','w')
	file1_dic = {}
	file2_dic = {}
	dmel_boleae = {'X':'I', '3R':'II','2R':'III','3L':'IV','2L':'Un','4':'Un','MT':'Un','Y':'Un','#Replicon':'Un'}
	a=b=c=d=e=f=g=z= 0
	for line in open(file1,'r'):
		sline = line.strip().split()
		geneID = sline[0]
		if geneID not in file1_dic.keys():
			file1_dic.setdefault(geneID,sline[1])
		elif geneID in file1_dic.keys():
			sys.stderr.write(line)
	print('\nMoving on...\n')
	oline = open(file2,'r')
	next(oline)
	for line in oline:
		sline = line.strip().split()
		try:
			geneID = sline[0].split("|")[3]
		except:
			geneID = sline[0].strip()
		locus = sline[13]
		if locus not in file2_dic.keys() and geneID in file1_dic.keys():
			file2_dic.setdefault(locus,'\t'.join([sline[0],sline[1],sline[13],sline[15],file1_dic[geneID]]))
			a += 1
		elif geneID not in file1_dic.keys():
			newfile2.write(line)
	print(list(file1_dic.items())[0])
	print(list(file2_dic.items())[0])
	print('\nNow querrying...\n')
	for line in open(file3,'r'):
		id1 = line.strip().split()[0]
		id2 = line.strip().split()[6]
		if id2 in file2_dic.keys():
			to_print = '\t'.join([file2_dic[id2],id1,dmel_boleae[id1]])+'\n'
			newfile.write(to_print)
			z += 1
		else:
			b += 1
			#sys.stderr.write(line)
	newfile.close()
	print('\n'+str(100*a/len(file1_dic.keys())))
	print('\n'+str(100*z/len(file1_dic.keys())))
	print('\n'+str(z)+'\t'+str(b))
	
def find_common_entries10(file1, file2,file3): #gtf, fusions, bed 1 is genome, 2 is TE names
	newfile = open(file3+'_edited','w')
	newfile2 = open(file3+'_Errors','w')
	newfile3 = open(file3+'_KeyErrors','w')
	#newfile.write('\t'.join(['contig','length','chromosome','No._genes\n']))
	#find_common_entries10('contig_assignmentbyprobes.sorted', 'gene_contigs.sorted' , '10x_MP_ONT_PacBio.0.5.pilon2_readlengths.txt')
	file1_dic = {}
	file2_dic = {}
	a=b=c=d=e=f=g=z= 0
	for line in open(file1,'r'):
		sline = line.strip().split()
		geneID = sline[0]
		if geneID not in file1_dic.keys():
			#file1_dic.setdefault(geneID,line.strip())
			file1_dic.setdefault(sline[2],geneID)
		elif geneID in file1_dic.keys():
			sys.stderr.write(line)
		a+=1
	print('\nMoving on...\n')
	for line in open(file2,'r'):
		sline = line.strip().split() #sline = line.strip().split(':')
		geneID = sline[0] #geneID = sline[0]
		Class = sline[-1] #Class = ':'.join(sline[1:])
		if geneID not in file2_dic.keys():
			file2_dic.setdefault(geneID,Class)
			#a += 1
		else:
			sys.stderr.write(line)
		b+=1
	print(list(file1_dic.items())[0])
	print(list(file2_dic.items())[0])
	print('\nNow querrying...\n')
	for line in open(file3,'r'):
		try:
			if line.startswith('##g'):
				newfile.write(line)
				#continue
			elif line.startswith('##s'):
				newfile.write(line.replace(line.strip().split()[1],file1_dic[line.strip().split()[1]]))
				#continue
			else:
				sline = line.strip().split('\t')
				geneID = sline[0]
				target = sline[-1].split()[0].split(';')[-1].replace('Target=','')
				#print(geneID + '\t' + target)
				if geneID in file1_dic.keys() and target in file2_dic.keys():
					#to_print = '\t'.join([file1_dic[id1],file2_dic[id1]])+'\n'
					#to_print = '\t'.join([id1,leng,file1_dic[id1],file2_dic[id1]])+'\n' #contig assignment by blast and mapping
					newGeneID = file1_dic[geneID]
					newTarget = file2_dic[target].split(':')[0]
					newTarget2 = file2_dic[target] #file2_dic[target].split(':')[0]
					#Class = ':'.join(file2_dic[target].split(':')[1:])
					newfile.write(line.strip().replace(geneID,newGeneID).replace('Target='+target,'Target='+newTarget2).replace(target,newTarget)+'\n')
					#newfile.write(line.replace(geneID,newGeneID).replace(target,newTarget)) #newfile.write(line.strip().replace(geneID,newGeneID).replace(target,newTarget)+';'+ Class+'\n')
					c += 1
				#elif id1 not in file1_dic.keys() and id1 in file2_dic.keys():
				#	to_print = '\t'.join([id1,leng,'Un',file2_dic[id1]])+'\n'
				#	newfile.write(to_print)
				#	c += 1
				else:
					#to_print = '\t'.join([id1,leng,'Un','0'])+'\n'
					newfile2.write(line)
					d += 1
					#sys.stderr.write(line)
				e+=1
		except:
			newfile3.write(line)
			f+=1
	newfile.close()
	newfile2.close()
	newfile3.close()
	print('\n'+str(100*c/len(file1_dic.keys())))
	print('\n'+str(100*c/len(file2_dic.keys())))
	print('\n'+str(100*c/e))
	print('\n'+'\t'.join([str(a),str(b),str(c),str(d),str(e),str(f)]))
	
def find_common_entries11(file1, file2):
	newfile = open(file2+'.bed','w')	
	file1_dic = {}
	w=x=y=z = 0
	for line in open(file1,'r'):
		sline = line.strip().split('\t')
		key = sline[0].split("_")[0].replace('mir','miR')
		if key not in file1_dic.keys():
			#file1_dic.setdefault(key, [sline[0],sline[1],sline[8],sline[9]])
			file1_dic.setdefault(key, [sline[0],sline[1],str(min(int(sline[8]),int(sline[9]))), str(max(int(sline[8]),int(sline[9])))])
		elif geneID in file1_dic.keys():
			sys.stderr.write(line)
	print('\nMoving on...\n')
	print('Now querrying...')
	for line in open(file2,'r'):
		if line.startswith('read_name'):
			continue
		sline = line.strip().split('\t')
		key = sline[0].split("_")[0]
		ref = sline[1]
		coord1 = str(min(int(sline[26]), int(sline[27])))
		coord2 = str(max(int(sline[26]), int(sline[27])))
		if key.replace('-3p','').replace('-5p','') in file1_dic.keys():
			items = file1_dic[key.replace('-3p','').replace('-5p','')]  #items = file1_dic[key.rstrip('-3p5p')] 
			#if int(str(int(items[2]) - int(coord1)).strip('-')) <= 100 and float(sline[8]) <= 1 and items[1] == sline[1]:
			if int(items[2]) <= int(coord1) and (int(items[3])+10) >= int(coord2) and float(sline[8]) <= 1 and items[1] == sline[1]:
				to_print = '\t'.join([sline[0],sline[1],coord1,coord2,items[0],items[2],items[3]]) + '\n'
				newfile.write(to_print)
				w+=1 #mature miRNA aligned to genome and captured well in blastn
			else:
				#print('\t'.join([sline[0],sline[1],coord1,coord2,'NN',items[0],items[1],items[2],items[3]]))
				x += 1 #mature miRNA aligned to genome and also in blastn but not well captured well in blastn
			z += 1
		else:
			y+=1 #mature miRNA aligned to genome but not in blast result
	print('\n'+str(100*z/len(file1_dic.keys())))
	print('\n'+str(z))
	print([w,x,y,z])

def find_common_entries12(file1, file2,file3): #gtf, bed, fusions,  1 is genome, 2 is TE names
	newfilex = open('readsNotStrandedExtra','w')
	#newfile1 = open(file3+'_readsNotStranded','w')
	#newfile2 = open(file3+'_singleExonNonStrandedOppDir','w')
	#newfile3 = open(file3+'_singleExonNonStrandedSameDir','w')
	#newfile4 = open(file3+'_SameDir','w')
	file1_dic = {}
	file2_dic = {}
	a=b=c=d=e=f=g=h=z= 0
	for line in open(file1,'r'):
		sline = line.strip().split() #[-1].split(';')
		if len(sline) <3 or line.startswith('#'):
			continue
		else:
			if sline[2] == 'exon':
				gene_IDPat = re.compile(r'gene_id "(.*?)"')
				gene_ID = gene_IDPat.search(line).group(1)
				transcpt_idPat = re.compile(r'transcript_id "(.*?)"')
				transcpt_id = transcpt_idPat.search(line).group(1)
				
				if gene_ID not in file1_dic.keys():
					file1_dic.setdefault(gene_ID,[transcpt_id])
				else:
					file1_dic[gene_ID].append(transcpt_id)
		a+=1
	print('\nMoving on...\n')
	for line in open(file2,'r'):
		if line.split()[3] not in file2_dic.keys():
			file2_dic.setdefault(line.split()[3],[line.split()[1]]+[line.split()[5]]+line.strip().split()[-3:])
		b+=1
	print(list(file1_dic.items())[0])
	print(list(file2_dic.items())[0])
	print('\nNow querrying...\n')
	for line in open(file3,'r'):
		try:
			gene1 = line.split()[2].split('--')[0]
			gene2 = line.split()[2].split('--')[1]
			transcpt1 = list(set(file1_dic[gene1]))
			transcpt2 = list(set(file1_dic[gene2]))
			#print(transcpt1)
			#print(transcpt2)
			for trans in range(len(transcpt1)):
				trans1_det = file2_dic[transcpt1[trans]]
				trans1_Strt = int(trans1_det[0])
				trans1_dir = trans1_det[1]
				trans1_exons = trans1_det[2]
				trans1_exoStrt = trans1_det[-1].strip('\n,').split(',')
				#print(trans1_exoStrt)
				for trans2 in range(len(transcpt2)):				
					trans2_det = file2_dic[transcpt2[trans2]]
					trans2_Strt = int(trans2_det[0])
					trans2_dir = trans2_det[1]
					trans2_exons = trans2_det[2]
					trans2_exoStrt = trans2_det[-1].strip('\n,').split(',')
					if trans1_dir != trans2_dir and int(trans1_exons) > 1 and trans1_exons == trans2_exons:
						exon_equal = 0
						for i in range(len(trans1_exoStrt)-1):
							i =+ 1
							if (int(trans2_exoStrt[i])+trans2_Strt)-10 <= (int(trans1_exoStrt[i])+trans1_Strt) <= (int(trans2_exoStrt[i])+trans2_Strt)+10:
								exon_equal += 1
						if exon_equal+1 == int(trans1_exons):
							newfile1.write(line)
							d += 1 #Delete all novel genes in this category where the novel gene is fused with NCBI gene
					elif trans1_dir != trans2_dir and int(trans1_exons) > 1 and trans1_exons != trans2_exons:
						for i in range(len(trans1_exoStrt)-1):
							i =+ 1
							if (int(trans2_exoStrt[i])+trans2_Strt)-10 <= (int(trans1_exoStrt[i])+trans1_Strt) <= (int(trans2_exoStrt[i])+trans2_Strt)+10:
								newfilex.write(line)
								continue
					elif trans1_dir != trans2_dir and int(trans1_exons) == 1 and trans1_exons == trans2_exons and line.split()[3] == 'Yes':
						newfile2.write(line)
						e += 1 #Remove the ones where NCBI gene is fused with novel gene. For fused novel gene retain only
					elif trans1_dir == trans2_dir and int(trans1_exons) == 1 and trans1_exons == trans2_exons and line.split()[3] == 'Yes':
						newfile3.write(line)
						f += 1
						#Delete all novel genes in this category where the novel gene is fused with NCBI gene
						#For novel genes that overlap each other, get expression data and remove the lower expressed one
					elif trans1_dir == trans2_dir and line.split()[3] == 'Yes':
						newfile4.write(line)
						g += 1

		except:
			#newfile3.write(line)
			h+=1
	#newfile1.close()
	#newfile2.close()
	newfilex.close()
	print(d)
	print(e)
	print(f)
	print(g)
	#print('\n'+str(100*c/len(file1_dic.keys())))
	#print('\n'+str(100*c/len(file2_dic.keys())))
	#print('\n'+str(100*c/e))
	#print('\n'+'\t'.join([str(a),str(b),str(c),str(d),str(e),str(f)]))

def find_common_entries13(file1, file2):
	newfile = open(file2+'.good.gtf','w')
	newfileBad = open(file2+'.bad.gtf','w')
	file1_list = []
	a=b=c=d=e = 0
	for line in open(file2,'r'):
		file1_list += [line.strip()]
		
	print('\nMoving on...\n')
	for line in open(file1,'r'):
		gene_IDPat = re.compile(r'gene_id "(.*?)"')
		gene_ID = gene_IDPat.search(line).group(1)

		if not gene_ID in file1_list:
			newfile.write(line)
			a += 1
		else:
			newfileBad.write(line)
			b += 1
	newfile.close()
	newfileBad.close()
	print(str(a)+'\t'+str(b))
	
def find_common_entries14(file1, file2):
	newfile = open(file1+'.10x_74X_Y-contigs','w')	
	file1_dic = {}
	missing = []
	a=b=c=w=x=y=z = 0
	for line in open(file1,'r'):
		sline = line.split('\t')
		key = sline[0]
		if key not in file1_dic.keys():
			file1_dic.setdefault(key, line)
	print('\nMoving on...\n')
	print('Now querrying...')
	for line in open(file2,'r'):
		if line.startswith('>'):
			key = line.strip('>\n')
			a += 1
			if key in file1_dic.keys():
				b += 1
				newfile.write(file1_dic[key])
			else:
				missing += [key]
	print([a,b])
	print(missing)
	
def find_common_lines(file1, file2):
	newfile = open(file1+'_common_positive', 'w')
	newfile2 = open(file1+'_negative', 'w')
	ofile1 = open(file1)
	lfile1 = ofile1.readlines()
	ofile2 = open(file2)
	lfile2 = ofile2.readlines()
	k = 0
	n = 0
	m = 0
	bar = progressbar.ProgressBar(maxval=len(lfile2), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in lfile2:
		if line.split('\t')[0]+'\n' in lfile1:
			n += 1			
			newfile.write(line)
			lfile1.remove(line.split('\t')[0]+'\n')
		else:
			newfile2.write(line)
			k += 1
		m += 1		
		bar.update(m)
	bar.finish()
	newfile.close()
	newfile2.close()
	#print('There are ' + str(len(lfile1)) + ' lines in ' + file1 +
	print(str(n) + '\t' + str(k))
	
def find_common_lines2(file1, file2):
	newfile = open(os.path.dirname(file1)+'/'+os.path.basename(file1)+'.resorted.bed', 'w')
	ofile1 = open(file1)
	lfile1 = ofile1.readlines()
	ofile2 = open(file2)
	lfile2 = ofile2.readlines()
	k=m=n=o=p=q=v=w = 0
	print([lfile1[0].strip()])
	print([lfile2[0].split()[3]])
	bar = progressbar.ProgressBar(maxval=len(lfile2), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line1 in lfile1:
		for line2 in lfile2:
			if line1.strip() == line2.split()[3]:
				newfile.write(line2)
				k += 1
			elif line1.strip() == line2.strip():
				n += 1
				p += 1
			elif line1.strip().strip("-like") in line2.strip():
				o += 1
				p += 1
				v += 1
				
		if p == 0:
			q += 1
			#newfile.write(line1)
		else:
			p = 0
		if v > 0:
			w += 1
			v == 0
		m += 1		
		bar.update(m)
	bar.finish()
	newfile.close()

	print('Exact matching proteins are ' + str(n))
	print('Non exact matching proteins are ' + str(w))
	print('Total hits for non exact matching proteins are ' + str(o))
	print('Missing proteins are ' + str(q))
	print('\nDiscovered ' + str(k) + ' genes' + ' which is %s percent of all the names u gave' %(str("{0:.1f}".format(round(100*k/len(lfile1)),1))))

def remove_lines(file,fasta):
	dirc = os.path.dirname(fasta)
	basename = os.path.basename(fasta)
	newfile = open(dirc+'/'+basename.replace('fasta','2.fasta'),'w')
	#newfile2 = open(dirc+'/'+basename.replace('.fasta','_repeatlines.fasta'),'w')
	
	a=b=c=d=e=f=g=h=j=k=l=n=m=x = 0
	
	file1 = open(file,'r')
	file1 = file1.readlines()
	
	file2 = open(fasta, 'r')
	length = int(os.popen('wc -l %s' %(fasta)).read().split()[0])
	print('We are looping over ' + str(length/2) + ' fasta sequences')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x<length:
		#'''
		name=file2.readline()
		seq =file2.readline()	
		if name.strip(">\n").split()[0]+'\n' not in file1:
		#if name not in file1:
			newfile.write(name+seq)
			a += 1
		else:
			#newfile2.write(name+seq)
			b += 1
		x+=2
		'''
		name=file2.readline()
		if name.strip(">\n").split()[0]+'\n' not in file1:
			newfile.write(name)
			a += 1
		else:
			b += 1
		x += 1
		'''
		bar.update(x)
	bar.finish()
	newfile.close()
	#newfile2.close()
	print('Found ' + str(a) + ' unique fasta sequences, and ' + str(b) + ' repeated fasta sequences')
	
def find_and_manip_matching_lines(file1, file2): #file2 is the larger one, we store it in memory and query it as we loop over file1
	###creat dictionary of file2
	dirc = os.path.dirname(file1)
	file2_dict = {}	
	
	x = 0
	
	length = int(os.popen('wc -l %s' %(file2)).read().split()[0]) #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
	length2 = length
	ofile2=open(file2,'r')
	#nfilex = open(dirc+'/'+ 'shit', 'w')
	next(ofile2)
	print('\nCreating dictionary for %s lines\n' %(str(length)))
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length: #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
		'''
		a=ofile2.readline()
		b=ofile2.readline()
		d=a.split("|")[0]
		e=a.split("|")[1]		
		file2_dict.setdefault(d.strip('>')+'_'+e, [a,b])		
		x+=2
		bar.update(x)
		'''
		'''
		a=ofile2.readline()
		#if a.startswith('PB'):
		#	b=a.strip().split('_')[0].strip(".path1")
		b=a.strip().split('_')[0].replace(".path1","")
		file2_dict.setdefault(b,a)
		x+=1
		'''
		''''
		a=ofile2.readline()
		a=a.split()
		if len(a) >=5:
			b=a[-1].strip() #.split()[4].strip(">") #b=a.strip().split()[2].strip(">") #
			file2_dict.setdefault(b,'\t'.join([a[0],a[1],b])) #file2_dict.setdefault(b,a.split()[1]) # #
		x+=1
		bar.update(x)
		''' #This one is to get homologues
		a=ofile2.readline()
		b=a.split()
		try:		
			#c=b[0].split("|")[0]
			d=b[1].split("|")[1]
			if not d in file2_dict.keys():
				file2_dict.setdefault(d,[a])
			else:
				file2_dict[d].append(a)
		except IndexError:
			print(a)
			#nfilex.write('\t'.join(b))
			#print('IndexError: list index out of range')
		x+=1
		bar.update(x)
	bar.finish()
	#nfilex.close()
	print(list(file2_dict.items())[0])
	print('\nDone creating dictionary, moving on to querrying...\n')	
	
	dirc = os.path.dirname(file1)
	name = os.path.basename(file2)
	nfile1 = open(dirc+'/'+ name + '.plus.dmel.names', 'w') #open(name.replace('_2','.gbk.fasta'), 'w')
	#nfile1.write(file2_dict['X1H'])
	nfile2 = open(dirc+'/'+ name + '.dmelnames.not.identified', 'w') #open(name.replace('_2','.gbk.nofasta'), 'w')
	
	a=b=c=d=e=f=g = 0
	
	length = int(os.popen('wc -l %s' %(file1)).read().split()[0])
	print('We are querying dictionary for ' + str(length) + ' items\n')
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for line in open(file1):
		'''
		gene_id = line.strip().split()[2]+'_'+line.strip().split()[1]
		if b < 1:
			print([gene_id])
		if gene_id in file2_dict.keys():
			seq_ID = file2_dict[gene_id][0]
			seq_ID = seq_ID.replace('>','>'+line.strip().split()[0]+'|')
			nfile1.write(seq_ID+file2_dict[gene_id][1])
			#nfile1.write('\t'.join([line.strip(),file2_dict[gene_id]]))
			a += 1
		else:
			nfile2.write(line)
		b+=1		
		bar.update(b)
		'''
		'''
		if line.strip().split()[1] in file2_dict.keys():
			nfile1.write('\t'.join([line.strip().split()[1],file2_dict[line.strip().split()[1]]]))
			#del file2_dict[line.strip()]		
			a += 1
		else:
			nfile2.write(line)
			b += 1
		c+=1		
		bar.update(c)
		'''
		'''
		if line.strip() in file2_dict.keys():
			item = file2_dict[line.strip()].split()
			filter(lambda ax: ax != '', item)
			nfile1.write('\t'.join(item)+'\n') #nfile1.write(file2_dict[line.strip()]) #			
			del file2_dict[line.strip()]		
			a += 1
		else:
			nfile2.write(line)
			b += 1
		'''
		'''
		gene_id = line.split()[1].split('|')[1]
		if gene_id in file2_dict.keys():
			nfile1.write('\t'.join([line.strip(),file2_dict[gene_id]])+'\n')
			a += 1
		c+=1		
		bar.update(c)
		'''
		line2 = line.split('\t')[-1]
		if line2.strip() in file2_dict.keys():
			for i in file2_dict[line2.strip()]:
				i = '\t'.join(i.split())
				nfile1.write(i.strip() +'\t'+ '\t'.join(line.split('\t')[:3]) +'\t'+ line.split('\t')[-1])
			#nfile1.write(line.split()[0]+'\t'+file2_dict[line.split()[0]]+'\n')
				a += 1
			del file2_dict[line2.strip()]
		else:
			nfile2.write(line)
			b += 1
		c+=1
		bar.update(c)
	bar.finish()
	print(file2_dict.keys())
	'''
	for x,y in file2_dict.items():
		nfile1.write(y)
	'''
	nfile1.close()
	nfile2.close()
	print('\n'+str(a)+'\t'+str(b))
	print('\nFound ' + str(a) + ' matches which is %s percent of all the queries u gave...\n' %(str("{0:.1f}".format(round(100*a/length2),1))))
	#print('\nFound ' + str(a) + ' matches which is %s percent of all the queries u gave...\n' %(str("{0:.1f}".format(round(100*a/len(file2_dict.keys())),1))))

def find_and_manip_matching_lines2(file1):
	dics = []
	x = 0

	dirc = os.path.dirname(file1)
	for i in files: #range(len(files)):
		dictn = i+'_dict'
		dictn = {}
		dics += [dictn]
		file_x = dirc+'/'+ files[i]
		length = int(os.popen('wc -l %s' %(file_x)).read().split()[0]) #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
		file_x=open(file_x,'r')
		print('\nCreating dictionary for %s lines\n' %(str(length)))
		bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
		bar.start()
		while x < length:
			a=file_x.readline()
			b=a.split()		
			dictn.setdefault(b[3],'\t'.joint([b[5],b[6].strip()]))
			x+=1
			bar.update(x)
		bar.finish()
	print('\nDone creating dictionary, moving on to querrying...\n')	
	#'''
	r = 0
	#print(len(file2_dict.keys()))
	for x,y in file2_dict.items():
		if r < 1:
			print([x])#,y)
			r +=1
	#'''
	
	dirc = os.path.dirname(file1)
	name = os.path.basename(file1)
	nfile1 = open(dirc+'/'+ name + '.Bo_homologue', 'w') #open(name.replace('_2','.gbk.fasta'), 'w')
	nfile1.write('\t'.join(['term.name','MD_adjust_p_value','MD_overlap_perc','MS_adjust_p_value','MS_overlap_perc',
							'TOP_adjust_p_value','TOP_overlap_perc','ZYG_adjust_p_value','ZYG_overlap_perc'])+'\n')
	
	a=b=c=d=e=f=g = 0
	MD_adjust_p_value ='0'
	MD_overlap_perc='0'
	MS_adjust_p_value='0'
	MS_overlap_perc='0'
	TOP_adjust_p_value='0'
	TOP_overlap_perc='0'
	ZYG_adjust_p_value='0'
	ZYG_overlap_perc='0'
	for i in dics:
		other_dics = dics.remove('i')
		dictn = file+str(i)+'_dict'		
		for x,y in i.items(): #for j in other_dics:
			term_name = x
			if 'MD' in i:
				MD_adjust_p_value = y[0]
				MD_overlap_perc=y[1]
			elif 'MS' in i:
				MS_adjust_p_value=y[0]
				MS_overlap_perc=y[1]
			elif 'TOP' in i:
				TOP_adjust_p_value=y[0]
				TOP_overlap_perc=y[1]
			elif 'ZYG' in i:
				ZYG_adjust_p_value=y[0]
				ZYG_overlap_perc=y[1]
			#for j in other_dics: #for x,y in i.items():
				#if x in j.keys():
					#if 'MD' in j:
					#	MD_adjust_p_value = j[x].[0]
					#	MD_overlap_perc=j[x].[1]
					#elif 'MS' in j:
					#	MS_adjust_p_value=j[x].[0]
					#	MS_overlap_perc=j[x].[1]
					#elif 'TOP' in j:
					#	TOP_adjust_p_value=j[x].[0]
					#	TOP_overlap_perc=j[x].[1]
					#elif 'ZYG' in j:
					#	ZYG_adjust_p_value=j[x].[0]
					#	ZYG_overlap_perc=j[x].[1]
			#nfile1.write('\t'.join([term_name,MD_adjust_p_value,MD_overlap_perc,MS_adjust_p_value,MS_overlap_perc,
			#						TOP_adjust_p_value,TOP_overlap_perc,ZYG_adjust_p_value,ZYG_overlap_perc])+'\n')
	#nfile1.close()	
def remove_funkylines(fasta):
	dirc = os.path.dirname(fasta)
	basename = os.path.basename(fasta)
	print(basename)
	try:
		fasta = check_fasta.check_fasta_fmt(fasta)
	except:
		if not os.path.exists(fasta):
			sys.exit(-1)
	newfile = open(dirc+'/'+basename.replace('fasta','edited.fasta'),'w')
	length = int(os.popen('wc -l %s' %(fasta)).read().split()[0]) #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
	file_x=open(fasta,'r')
	c=x=k=0
	empty = []
	bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length:
		a=file_x.readline()
		b=file_x.readline()
		if len(b.strip().strip("N")) < 100:
			empty+=[a]
		elif len(b.strip().strip("N")) >= 100:
			newfile.write(a+b.strip().strip("N").replace('te','').strip('e')+'\n')
		x+=2
		c+=1
		bar.update(c)
	bar.finish()
	newfile.close()
	print(empty)
	
def get_dmel_boleae_orthos(file1, file2):
	dirc = os.path.dirname(file1)
	name = os.path.basename(file1)
	newfile = open(dirc+'/'+name + '.dmelOrthos','w')

	file2_dic = {}
	z = 0

	#print(list(file1_dic.items())[0])
	print('\nMoving on...\n')	
	for line in open(file2,'r'):
		sline = line.strip().split()
		geneID = sline[0].split("|")[0]
		dmel = sline[1].split("|")[1]
		if geneID not in file2_dic.keys():
			file2_dic.setdefault(geneID,dmel)
		elif geneID in file2_dic.keys():
			sys.stderr.write(line)
	#print(list(file2_dic.items())[0])

	for y in open(file1, 'r'):
		if y.split()[1] in file2_dic.keys():
			newfile.write(file2_dic[y.split()[1]]+'\n')
			z += 1
		else:
			print(y)
			#newfile.write(y)
	
	print('\n'+str(100*z/len(file2_dic.keys())))
	
def find_gff3_genes(gff3, file2):
	dirc = os.path.dirname(file2)
	nfile1 = open(dirc+'/'+'zygotic_early_genes_coordinates', 'w')
	
	c=d=e=f=g= 0
	genes = {}
	length = int(os.popen('wc -l %s' %(gff3)).read().split()[0])
	print('We are looping over ' + str(length) + ' lines')
	file1=open(gff3,'r')
	lfile1=file1.readlines()
	bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for i in lfile1:
		ftrs = i.split('\t')
		if len(ftrs) > 3 and ftrs[2]=='gene':
			GeneID = "."
			description = "."
			gene_length = int(ftrs[4]) - int(ftrs[3]) +1
			ftrs2 = ftrs[-1].split(";")
			gene_name=ftrs2[0].strip("ID=")
			if gene_name not in genes.keys():
				genes.setdefault(gene_name,'\t'.join([gene_name,ftrs[0]+':'+ftrs[3]+'-'+ftrs[4],str(gene_length)+'\n']))
			else:
				print(gene_name)		
		d+=1
		bar.update(d)
	bar.finish()
	#print(list(genes.items())[0])
	#file1=open(file2,'r')
	#lfile1=file1.readlines()
	#print(lfile1[10].split('\t')[0])
	for ix in open(file2, 'r'):
		ftrs = ix.split('\t')
		if ftrs[0] in genes.keys():
			nfile1.write(genes[ftrs[0]])
			g += 1
	nfile1.close()
	print(g)
	
def find_longest_transcript_gtf(gtf):
	dirc = os.path.dirname(gtf)
	name = os.path.basename(gtf)
	newfile = open(dirc+'/'+name + '.longest_transcript','w')
	newfile.write('\t'.join([ 'gene_id','longest_transcript\n' ]))
	gene_dic = {}
	transcpt_dic = {}
	a=b = 0
	
	print('\nMoving on...\n')
	for line in open(gtf,'r'):
		sline = line.strip().split() #[-1].split(';')
		if len(sline) <3 or line.startswith('#'):
			continue
		else:
			if sline[2] == 'exon':
				exon_len = int(sline[4]) - int(sline[3])
				gene_IDPat = re.compile(r'gene_id "(.*?)"')
				gene_ID = gene_IDPat.search(line).group(1)
				transcpt_idPat = re.compile(r'transcript_id "(.*?)"')
				transcpt_id = transcpt_idPat.search(line).group(1)
				
				if gene_ID not in gene_dic.keys():
					gene_dic.setdefault(gene_ID,[transcpt_id])
				else:
					gene_dic[gene_ID].append(transcpt_id)
				if transcpt_id not in transcpt_dic.keys():
					transcpt_dic.setdefault(transcpt_id,[exon_len])
				else:
					transcpt_dic[transcpt_id].append(exon_len)
	for gene, transList in gene_dic.items():
		longest_trans = ''
		trans_len = 0
		for i in list(dict.fromkeys(transList)):
			transx_len = sum(transcpt_dic[i])
			if transx_len > trans_len:
				longest_trans = i
		newfile.write(gene+'\t'+longest_trans+'\n')
	newfile.close()
	
def gft_fusion(fusions,gtf):
	dirc = os.path.dirname(fusions)
	name = os.path.basename(fusions)
	newfile = open(dirc+'/'+name + '.new','w')
	newfile.write('\t'.join([ 'read','No.Fused','fusions','IfFused','Fusion1_ref','Fusion2_ref\n' ]))
	gene_dic = {}
	for line in open(gtf,'r'):
		sline = line.strip().split()
		if len(sline) <3 or line.startswith('#'):
			continue
		else:
			if sline[2] == 'exon':
				ref = sline[0]
				gene_IDPat = re.compile(r'gene_id "(.*?)"')
				gene_ID = gene_IDPat.search(line).group(1)
				
				if gene_ID not in gene_dic.keys():
					gene_dic.setdefault(gene_ID,ref)
				else:
					yme = 0
	for line in open(fusions,'r'):
		sameRef = ''
		sline = line.strip().split()
		fusdGenes = sline[-1].split('--')
		if gene_dic[fusdGenes[0]] == gene_dic[fusdGenes[1]]:
			sameRef = 'Yes'
		else:
			sameRef = 'No'
		newfile.write('\t'.join([ line.strip(),sameRef,gene_dic[fusdGenes[0]],gene_dic[fusdGenes[1]]+'\n' ]))
	newfile.close()
	
def remove_bad_novel_genes(fusions, tpm):
	dirc = os.path.dirname(fusions)
	name = os.path.basename(fusions)
	newfile1 = open('NovelGenesToKeep','w')
	tpm_dic = {}
	a=b = 0

	for line in open(tpm,'r'):
		sline = line.strip().split()
		if sline[0].split('_')[0] not in tpm_dic.keys():
			tpm_dic.setdefault(sline[0].split('_')[0],float(sline[1]))
		elif sline[0] in tpm_dic.keys():
			sys.stderr.write(line)			

	print(list(tpm_dic.items())[0])
	print('\nMoving on...\n')
	for line in open(fusions,'r'):
		del tpm_dic[line.strip()]
	for novelgene in tpm_dic.keys():
		newfile1.write(novelgene+'\n')
	'''
	for line in open(fusions,'r'):
		sline = line.strip().split('--')
		gene1_tpm = tpm_dic[sline[0]]
		gene2_tpm = tpm_dic[sline[1]]
		if gene1_tpm < gene2_tpm:
			newfile1.write(sline[0]+'\n')
		else:
			newfile1.write(sline[1]+'\n')
	'''
	newfile1.close()
	
def find_some_quickMatches(file1, file2): #to remove, nvlGene
	#This script takes a fasta file and a file with the names of reads u want to exclude from the fasta. It will create a dictionary for the fasta seqs and a list for bad reads. 
	#If a name in the file is found in the fasta as well, the corresponding sequence and the name are not returned. Pay attention to the \n and > and fasta format, it has to be single line for seqs
	a=b=c=d=0
	print("\nCreating dictionary of the fasta...\n")
	lopenfile = []
	for line in open(file1,'r'):
		lopenfile += ['.'.join(line.split('.')[:2])]
	l = len(lopenfile)

	path = os.path.dirname(file2)
	name = os.path.basename(file2)
	outa=open(path+'/'+name+'.NOJUNK','w')
	outb=open(path+'/'+name+'.JUNK','w')
	
	for line in open(file2,'r'):
		line1 = '.'.join(line.split()[0].split('.')[:2])
		if line1 not in lopenfile:
			outa.write(line)
			a += 1
		else:
			outb.write(line)
		b+=1
	outa.close()
	print('Total reads = ' + str(b) + '\nPrinted reads = ' + str(a) + '\nBad reads = ' + str(l) + '\nUn-printed reads = '+str(b-a))
	#print(lopenfile2)	
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_5H/combined/assembly/tofu'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/ncbi_genome/dmel'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/ncbi_genome/dmel/de_renzis'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/ncbi_genome/full_annotation/proteins/functional_annotation'
dr = '/home/banthony/scratch/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/rna_velocity/hamed/annotation'
dr = '/home/banthony/scratch/analysis/illumina/olivefly/rna_seq/C010_11/novaseq_A00266_0083/Bo_All'
dr = '/home/banthony/scratch/analysis/nanopore/olivefly/ncbi_genome'
dr = '/home/banthony/scratch/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/SQANTI_tofu_min2_polyAtrim/sqanti_rpt'
dr1 = '/home/banthony/scratch/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/SQANTI_tofu_min2_polyAtrim/novel_genes/minimap2_find_wrong_novel_genes'
dr = '/home/banthony/scratch/analysis/nanopore/olivefly/ncbi_genome/full_annotation/proteins/functional_annotation'
dr = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/unmapped/sqanti'
dr2 = '/home/banthony/scratch/analysis/illumina/capitata/sra/rsem'

#find_common_entries(dr+'/'+'z.maternal_degraded_genes.Bo_homologue', dr+'/'+'exons_mandalorion_editted')
file1='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair_expn/C1-15H_counts_matrix_tpm_edited_bulk.tsv.averages'
file2='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/tofu/reflair/master_file'
file2='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair_expn/C1-15H_counts_matrix_bulk.tsv.tpm.tsv'
#find_common_entries2(file1,file2)
#find_common_entries2(dr+'/'+'z.maternal_degraded_genes.Bo_homologue', dr+'/'+'introns_mandalorion_editted')
#find_common_entries3(dr2+'/'+'ncbiR100.JAMgOGS1_minimap2_with_alignments.ed2.txt', dr+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds.gff3')
file1='/home/banthony/scratch/analysis/illumina/capitata/sra/rsem/C1H_pass_FL_canu.cutadapt.flair.isoforms.results_ed'
file1='/home/banthony/scratch/analysis/covid19/B004/samples_to_analyse/sample_lineage'
file2='/home/banthony/scratch/analysis/covid19/B004/samples_to_analyse/V4testSamples_test.csv'
#find_common_entries3(dr2+'/'+'CAJHJT01_C1-15H_pass_FL_canu.cutadapt.unmapped.flair.isoforms_ed.isoforms.results_ed', dr+'/'+'CAJHJT01_C1-15H_pass_FL.unmapped.flair.sqanti_classification.txt')
find_common_entries3(file1, file2)
#find_common_entries4(dr+'/'+'file1',dr+'/'+'file2',dr+'/'+'b.oleae_c.cap_gene_homologues')
#find_and_manip_matching_lines(dr+'/'+'fbgn_NAseq_Uniprot_fb_2018_03',dr+'/'+'Bo.all_proteins_UP000000803_blastp_sorted.txt')
#find_and_manip_matching_lines(dr+'/'+'Bo.all_proteins_UP000000803_blastp_sorted.txt' ,dr+'/'+'early_n_maternal.edited')
#find_and_manip_matching_lines(dr+'/'+'only_maternal',dr+'/'+'fbgn_fbtr_fbpp_fb_2018_03.tsv')
#find_and_manip_matching_lines(dr+'/'+'nor_early_or_maternal',dr+'/'+'fbgn_fbtr_fbpp_fb_2018_03.tsv')

#find_and_manip_matching_lines(dr+'/'+'only_early',dr+'/'+'fbgn_NAseq_Uniprot_fb_2018_03')
#find_and_manip_matching_lines(dr+'/'+'early_n_maternal',dr+'/'+'fbgn_NAseq_Uniprot_fb_2018_03')
#find_and_manip_matching_lines(dr+'/'+'only_maternal',dr+'/'+'fbgn_NAseq_Uniprot_fb_2018_03')
#find_and_manip_matching_lines(dr+'/'+'nor_early_or_maternal',dr+'/'+'fbgn_NAseq_Uniprot_fb_2018_03')

#find_and_manip_matching_lines(dr+'/'+'maternal_degraded_zygotic_new_ids.uniprot.Bo_homologue.gene',dr+'/'+'combined_mandalorion_absolute_ed')
#find_and_manip_matching_lines(dr+'/'+'maternal_zygotic_new_ids.uniprot.Bo_homologue.gene',dr+'/'+'combined_mandalorion_absolute_ed')
#find_and_manip_matching_lines(dr+'/'+'non_zygotic_secondary_new_ids.uniprot.Bo_homologue.gene',dr+'/'+'combined_mandalorion_absolute_ed')
#find_and_manip_matching_lines(dr+'/'+'zygotic_early_new_ids.uniprot.Bo_homologue.gene',dr+'/'+'combined_mandalorion_absolute_ed')
#find_and_manip_matching_lines(dr+'/'+'zygotic_primary_new_ids.uniprot.Bo_homologue.gene',dr+'/'+'combined_mandalorion_absolute_ed')
#find_and_manip_matching_lines(dr+'/'+'zygotic_pure_new_ids.uniprot.Bo_homologue.gene',dr+'/'+'combined_mandalorion_absolute_ed')
#find_and_manip_matching_lines(dr+'/'+'zygotic_pure_tc_new_ids.uniprot.Bo_homologue.gene',dr+'/'+'combined_mandalorion_absolute_ed')
#find_and_manip_matching_lines(dr+'/'+'zygotic_secondary_new_ids.uniprot.Bo_homologue.gene',dr+'/'+'combined_mandalorion_absolute_ed')
#find_and_manip_matching_lines(dr+'/'+'z.maternal_degraded_genes',dr+'/'+'Bo.all_proteins_UP000000803_blastp_sorted.txt')

#dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_allH/rna_velocity/zygotic_early'
#find_common_entries5(dr2+'/'+'coverage.1H.bed',dr+'/'+'zygotic_early_genes',dr+'/'+'full_splice_match_genes')
#find_common_entries5(dr2+'/'+'coverage.2H.bed',dr+'/'+'zygotic_early_genes',dr+'/'+'full_splice_match_genes')
#find_common_entries5(dr2+'/'+'coverage.3H.bed',dr+'/'+'zygotic_early_genes',dr+'/'+'full_splice_match_genes')
#find_common_entries5(dr2+'/'+'coverage.4H.bed',dr+'/'+'zygotic_early_genes',dr+'/'+'full_splice_match_genes')
#find_common_entries5(dr2+'/'+'coverage.5H.bed',dr+'/'+'zygotic_early_genes',dr+'/'+'full_splice_match_genes')
#find_common_entries5(dr2+'/'+'coverage.6H.bed',dr+'/'+'zygotic_early_genes',dr+'/'+'full_splice_match_genes')
#find_common_entries5(dr2+'/'+'coverage.female.bed',dr+'/'+'zygotic_early_genes',dr+'/'+'full_splice_match_genes')
#find_common_entries5(dr2+'/'+'coverage.male.bed',dr+'/'+'zygotic_early_genes',dr+'/'+'full_splice_match_genes')

#find_common_entries5(dr2+'/'+'coverage.1H.bed',dr+'/'+'maternal_degraded_genes_stages',dr+'/'+'full_splice_match_genes',dr+'/'+'busco_diptera_ancestral_genes')
#find_common_entries5(dr2+'/'+'coverage.2H.bed',dr+'/'+'maternal_degraded_genes_stages',dr+'/'+'full_splice_match_genes',dr+'/'+'busco_diptera_ancestral_genes')
#find_common_entries5(dr2+'/'+'coverage.3H.bed',dr+'/'+'maternal_degraded_genes_stages',dr+'/'+'full_splice_match_genes',dr+'/'+'busco_diptera_ancestral_genes')
#find_common_entries5(dr2+'/'+'coverage.4H.bed',dr+'/'+'maternal_degraded_genes_stages',dr+'/'+'full_splice_match_genes',dr+'/'+'busco_diptera_ancestral_genes')
#find_common_entries5(dr2+'/'+'coverage.5H.bed',dr+'/'+'maternal_degraded_genes_stages',dr+'/'+'full_splice_match_genes',dr+'/'+'busco_diptera_ancestral_genes')
#find_common_entries5(dr2+'/'+'coverage.6H.bed',dr+'/'+'maternal_degraded_genes_stages',dr+'/'+'full_splice_match_genes',dr+'/'+'busco_diptera_ancestral_genes')
#find_common_entries5(dr2+'/'+'coverage.female.bed',dr+'/'+'maternal_degraded_genes_stages',dr+'/'+'full_splice_match_genes',dr+'/'+'busco_diptera_ancestral_genes')
#find_common_entries5(dr2+'/'+'coverage.male.bed',dr+'/'+'maternal_degraded_genes_stages',dr+'/'+'full_splice_match_genes',dr+'/'+'busco_diptera_ancestral_genes')

#find_common_entries6(dr+'/'+'ancestral.pub_og_id', dr+'/'+'odb9v1_OG2genes.tab2', dr+'/'+'odb9v1_genes.tab2', dr+'/'+'Bo.all_proteins_UP000000803_blastp_sorted.txt.plus.dmel.names2')
#find_common_entries6(dr+'/'+'ancestral.pub_og_id', dr+'/'+'odb9v1_OG2genes.tab', dr+'/'+'odb9v1_genes.tab', dr+'/'+'Bo.all_proteins_UP000000803_blastp_sorted.txt.plus.dmel.names')

file1='/home/banthony/scratch/analysis/illumina/capitata/sra/rsem/C1H_pass_FL_canu.cutadapt.flair.isoforms.results_ed'
file2='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/flair.firstpass.q.counts_ed'
file3='/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/tofu/Ccap_flair.isoforms.tofu.collapsed.group.txt'
#find_common_entries7(file1,file2,file3)
#find_common_entries7('/home/banthony/scratch/analysis/combined/10x_All/polished_50X/10x_MP_ONT_PacBio.0.5.pilon2_readlengths.txt','Tag_and_sequence_names2017','cleaned_alignments.ed')
#find_common_entries8(dr+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gffread.gtf',dr1+'/'+'NovelGenesToKeep')
#find_common_entries9('gene_contigs', 'Bo.all_proteins_UP000000803_blastp_sorted.txt.plus.dmel.names', 'ProteinTable47_204923.txt')
#find_common_entries10('contig_assignmentbyprobes.sorted', 'gene_contigs.sorted' , '10x_MP_ONT_PacBio.0.5.pilon2_readlengths.txt')
#find_common_entries10('boleae_contigs_uniquely_placed', 'gene_contigs.sorted' , '10x_MP_ONT_PacBio.0.5.pilon2_readlengths.txt')
#find_common_entries10('Galaxy9-10x_MP_ONT_PacBio.0.5.pilon2_sl.noContaminants.ed2.fasta.fasta.matches', 'boleae_pirate_TEs_ed2.fasta.matches', 'Boleae_MU_v2_TEannot.gff3.gff3')
#find_common_entries10('Galaxy9-10x_MP_ONT_PacBio.0.5.pilon2_sl.noContaminants.ed2.fasta.fasta.matches', 'boleae_pirate_TEs_ed2.fasta.matches', 'zGalaxy304-Run1-TEannot.gff3')
#find_common_entries10('Galaxy9-10x_MP_ONT_PacBio.0.5.pilon2_sl.noContaminants.ed2.fasta.fasta.matches', 'boleae_pirate_TEs_ed2.fasta.matches', 'test')
#remove_funkylines(dr+'/'+'uncategorised_repeated_TEs.fasta')
#remove_funkylines(dr+'/'+'potentially_autonomous_TEs.fasta')
#remove_funkylines(dr+'/'+'non_autonomous_TEs.fasta')
#find_common_entries11(dr+'/'+'1',dr+'/'+'12')
#find_common_entries11('/home/banthony/scratch/analysis/combined/10x_All/polished_50X/miRNA/blast/insect_pre-miRNA_blastx.ed.txt','/home/banthony/scratch/analysis/combined/10x_All/polished_50X/miRNA/bowtie/MU_boleae_V2_insects.mature_miRNA_with_alignments.txt')
#get_dmel_boleae_orthos(dr+'/'+'adult','/home/banthony/scratch/analysis/nanopore/olivefly/ncbi_genome/full_annotation/proteins/functional_annotation/Bo.all_proteins_UP000000803_blastp_sorted.txt.plus.dmel.names')
##find_gff3_genes(dr1+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds.gff3',dr2+'/'+'zygotic_early_genes')
#find_gff3_genes(dr3+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gff3',dr2+'/'+'zygotic_early_genes')
#find_longest_transcript_gtf(dr3+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gffread.gtf')
#gft_fusion(dr3+'/'+'longest_transcript_6H.male.female_minimap2.fusions',dr2+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gffread.gtf')
#find_common_entries12(dr2+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gffread.gtf',
#					  dr2+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gff3_1.bed12',
#					  dr3+'/'+'longest_transcript_6H.male.female_minimap2.fusions.new.uniq') #gtf, bed, fusions 
#find_common_entries13(dr+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gffread.gtf',dr2+'/'+'final_list_of_novel_genes_to_remove')
#find_common_entries13(dr+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gff3_1.gencode.gtf',dr2+'/'+'final_list_of_novel_genes_to_remove')
#find_common_entries13(dr+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gffread.gtf',dr2+'/'+'final_list_of_novel_genes_to_remove')
#dr2+'/'+'ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gff3_1.bed12',
#dr2+'/'+'ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_scaffolds_plus_novel_genes_plus_10X_genes.gffread.gtf'
#dr3+'/'+'longest_transcript_6H.male.female_minimap2.fusions.new.uniq'
#remove_bad_novel_genes('final_list_of_novel_genes_to_remove','TPM_expressions')
#find_some_quickMatches(dr+'/'+'final_list_of_novel_genes_to_remove',dr+'/'+'novelGene')
dr1='/home/banthony/scratch/analysis/combined/10x_All/polished_50X/gnometognome'
dr='/home/banthony/scratch/analysis/combined/10x_All/y-chromosome'
#find_common_entries14(dr1+'/'+'LGAM02.1_olivefly_supernova_500kbc_74X_minimap2_with_alignments_edited.txt',dr+'/'+'olivefly_supernova_v1p1_corrvalid_group_gt10_500000bc_74X_seed0.pseudohap.1.repeatmasked.fillrepeat.Y.fasta')
dr1='/home/banthony/scratch/analysis/combined/10x_All/polished_50X/gnometognome/minimap2'
dr='/home/banthony/scratch/analysis/combined/10x_All/polished_50X/y_chromosome/ygs'
#find_common_entries14(dr1+'/'+'LGAM02.1_10x-ONT-PacBio_minimap2_with_alignments_edited.txt',dr+'/'+'YGS_validated_above80')

dr = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/Cap_all/flair/sqanti/newgtf'
file1='ERCC92_GCF_000347755.3_Ccap_2.1_genomic.gtf'
file2='C1-15H_pass_FL_canu.flair.sqanti_corrected.gtf.cds.gff'
file3='flair_FINAL_sqanti_plus_manually_filtered_isoforms'
#find_common_entries4(dr+'/'+file1, dr+'/'+file2, dr+'/'+file3)