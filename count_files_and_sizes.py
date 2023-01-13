#!/usr/bin/env python3
#we specifically need python >= 3.6 module load python/3.6.3
import re, os,sys #, check_fasta, humanize
import progressbar
from pathlib import Path
import find_fasta_seqs

def folders_in_path(path):
	##print('yea2')
	if not Path.is_dir(path):
		raise ValueError("argument is not a directory")
	yield from filter(Path.is_dir, path.iterdir())

def folders_in_depth(path, depth):
	##print('yea1')
	if 0 > depth:
		raise ValueError("depth smaller 0")
	if 0 == depth:
		yield from folders_in_path(path)
	else:
		for folder in folders_in_path(path):
			print(folder)
			yield from folders_in_depth(folder, depth-1)

def files_in_path(path):
	##print('yea3')
	if not Path.is_dir(path):
		raise ValueError("argument is not a directory")
	yield from filter(Path.is_file, path.iterdir())

def sum_file_size(filepaths):
	try:
		sum_x = sum([filep.stat().st_size for filep in filepaths])
	except:
		sum_x = 0
	return sum_x


	
#if __name__ == '__main__':
#	for folder in folders_in_depth(Path.cwd(),0):
#		##print('yes4')
#		#       vvvv quick hack to to use len(), does not perform well
#		files = list(files_in_path(folder))
#		total_size = sum_file_size(files)
#		print(f'{folder}: filecount:{len(files)}, total size:{humanize.naturalsize(total_size)}')

def sum_file_size(filepaths):
	return sum([os.stat(filep).st_size for filep in filepaths])

def count_files_and_sizes(rootDir):
	#newfile = open(rootDir+'/'+'files_and_folder_stats', 'w')
	newfile = open('/home/banthony/scratch/'+'files_and_folder_stats', 'w')
	for dirName, subdirList, fileList in os.walk(rootDir):
		files = []
		for name in fileList:
			file_x = os.path.abspath(os.path.join(dirName, name))
			if not os.path.islink(file_x):
				files += [file_x]
		total_size = sum_file_size(files)
		#newfile.write(f'{dirName}: subdirectories:{len(subdirList)}, filecount:{len(files)}, total_size:{humanize.naturalsize(total_size)}\n')
	newfile.close()

'''
def xcombine_fastq_files(rootDir):	
	for dirName, subdirList, fileList in os.walk(rootDir):
		files = []
		a = 0
		if len(subdirList) == 0 and len(fileList) > 0:
			for name in fileList:
				if name.endswith('.fq'):
					a+=1
					filex = name
			if len(fileList) == a:
				print(fileList)
				run_name = fileList[0].split('_')[1] #.replace('.fq','')
				cmd = 'cd %s && cat *.fq > %s' %(dirName,'runid_'+run_name+'.fastq')
				os.system(cmd)
				os.system('cd %s && rm *.fq' %dirName)
'''
def combine_fastq_files(rootDir):
	files = []
	file_paths = []
	for dirName, subdirList, fileList in os.walk(rootDir):		
		a = 0
		if len(subdirList) == 0 and len(fileList) >= 1:
			print('Now processing folder ' + dirName + ' Number of files = ' + str(len(fileList)))
			#run_name = fileList[0].split('_')[2]
			#path1 = os.path.abspath(os.path.join(dirName,'fastq_runid_'+run_name+'.fq'))
			#if os.path.exists(path1) and os.path.isfile(path1):
			#	os.unlink(path1)
			
			if len(list(Path(dirName).glob('*.fq'))) >= 1 and sum_file_size(list(Path(dirName).glob('*.fq'))) <= 100:
				for filepath in Path(dirName).glob('*.fq'):
					os.unlink(filepath)
			elif len(list(Path(dirName).glob('*.fq'))) >= 1 and len(list(Path(dirName).glob('*.fastq'))) == 0:
				number_of_reads = str( int(os.popen('wc -l %s' %( str(list(Path(dirName).glob('*.fq'))[0]) )).read().split()[0])/4 ).split('.')[0]
				files+=[dirName+'\t'+number_of_reads]
			elif len(list(Path(dirName).glob('*.fq'))) >= 1 and len(list(Path(dirName).glob('*.fastq'))) >= 1:
				#print(str(list(Path(dirName).glob('*.fq'))[0]).split('/')[-1].split('_')[2] + '\t' + str(list(Path(dirName).glob('*.fastq'))[0]).split('/')[-1].split('_')[2])
				if str(list(Path(dirName).glob('*.fq'))[0]).split('/')[-1].split('_')[2].strip('.fq') == str(list(Path(dirName).glob('*.fastq'))[0]).split('/')[-1].split('_')[2]:
					print('yes')
					if sum_file_size(list(Path(dirName).glob('*.fq'))) >= (0.9*sum_file_size(list(Path(dirName).glob('*.fastq')))):					
						for filepath in list(Path(dirName).glob('*.fastq')):
							os.unlink(filepath)
							#print('unlinked ' + str(list(filepath))
						number_of_reads = str( int(os.popen('wc -l %s' %( str(list(Path(dirName).glob('*.fq'))[0]) )).read().split()[0])/4 ).split('.')[0]
						files+=[dirName+'\t'+number_of_reads]
						
			if len(list(Path(dirName).glob('*.fq'))) == 0 and len(list(Path(dirName).glob('*.fastq'))) >= 1:
				file_paths = list(Path(dirName).glob('*.fastq'))
			#for name in fileList:
			#	if name.endswith('.fastq'):
			#		a+=1
			#		file_paths += [os.path.abspath(os.path.join(dirName, name))]
			if len(file_paths) > 1:
				run_name = fileList[0].split('_')[2]
				path1 = os.path.abspath(os.path.join(dirName,'fastq_runid_'+run_name+'.fq'))
				newfq = dirName+'/'+'fastq_runid_'+run_name+'.fq'
				#print(dirName+'/fastq_runid_'+run_name+'.fq')
				#cmd = 'cd %s && cat *fastq >> %s' %(dirName,'fastq_runid_'+run_name+'.fq')
				#os.system(cmd)
				#for name in fileList:
				#	cmd = 'cd %s && cat %s >> %s' %(dirName,name,'fastq_runid_'+run_name+'.fq')
				#	os.system(cmd)
				allreads = []
				for name in fileList:
					for eachline in open(dirName+'/'+name, 'r'):
						allreads += [eachline]
				newfile = open(newfq,'w')
				newfile.write(''.join(allreads))
				newfile.close()
				number_of_reads = str( int(os.popen('wc -l %s' %(newfq)).read().split()[0])/4 ).split('.')[0]
				files+=[dirName+'\t'+number_of_reads]
				#os.system('cd %s && rm *fastq' %dirName)
				if os.path.isfile(path1) and os.stat(path1).st_size >= (0.9*sum_file_size(file_paths)):
					#for filepath in file_paths:
					#	os.unlink(filepath)
					for filepath in Path(dirName).glob('*.fastq'):
						os.unlink(filepath)
					file_paths = []
	newfile = open(rootDir+'/basecalling_stats_2', 'w')
	newfile.write('\n'.join(files))
	newfile.close()

def rename_reads(rootDir):
	print('Renaming your reads now, fingers crossed...')
	directory = ''
	files = []
	for dirName, subdirList, fileList in os.walk(rootDir):		
		folder = dirName.split('/')[-1]
		if folder.startswith('barcode') or folder.startswith('unclassified'):
			print('Now processing folder: ' + dirName)
			if len(list(Path(dirName).glob('*.fq'))) == 0:
				continue
			elif len(list(Path(dirName).glob('*.fq'))) == 1 and len(list(Path(dirName).glob('*_renamed.fq'))) == 1:
				newfq = str(list(Path(dirName).glob('*.fq'))[0])
				number_of_reads = str( int(os.popen('wc -l %s' %(newfq)).read().split()[0])/4 ).split('.')[0]
				files += [newfq+'\t'+number_of_reads]
				continue
			elif len(list(Path(dirName).glob('*.fq'))) == 2:
				if str(list(Path(dirName).glob('*.fq'))[0]).split('/')[-1].split('_')[2].strip('.fq') == str(list(Path(dirName).glob('*.fq'))[1]).split('/')[-1].split('_')[2].strip('.fq'):
					if os.path.getsize(str(list(Path(dirName).glob('*.fq'))[0])) >= (0.9*os.path.getsize(str(list(Path(dirName).glob('*.fq'))[1]))):
						os.unlink(str(list(Path(dirName).glob('*.fq'))[0]).replace('_renamed',''))
						newfq = str(list(Path(dirName).glob('*.fq'))[0])
						number_of_reads = str( int(os.popen('wc -l %s' %(newfq)).read().split()[0])/4 ).split('.')[0]
						files += [newfq+'\t'+number_of_reads]
					else:
						print('Problematic folder = ' + dirName)
				else:
					print('Problematic folder = ' + dirName)
			elif len(list(Path(dirName).glob('*_renamed.fq'))) == 0 and len(list(Path(dirName).glob('*.fq'))) == 1:
				if folder.startswith('barcode'):
					directory = folder.replace('barcode','_BAR')
				elif folder.startswith('unclassified'):
					directory = ('_unclas')
				name = fileList[0] #for name in fileList:
				newfile = open(dirName+'/'+name.replace('.fq','_renamed.fq'), 'w')
				for line in open(dirName+'/'+name, 'r'):
					if line.startswith('@'):
						newfile.write(line.split()[0]+directory+'\t'+'\t'.join(line.split()[1:])+'\n')							
					else:
						newfile.write(line)
				newfile.close()
				newfq = dirName+'/'+name.replace('.fq','_renamed.fq')
				number_of_reads = str( int(os.popen('wc -l %s' %(newfq)).read().split()[0])/4 ).split('.')[0]
				files += [newfq+'\t'+number_of_reads]
				directory = ''				
				path1 = os.path.abspath(os.path.join(dirName,name.replace('.fq','_renamed.fq')))
				path2 = os.path.abspath(os.path.join(dirName,name))
				if os.path.isfile(path1) and os.stat(path1).st_size > (0.9*os.path.getsize(path2)):
					os.system('rm %s' %(path2))
		#elif folder.startswith('unclassified'):
		#	for name in fileList:
		#		if name.endswith('.fq'):
		#			files += [dirName+'/'+name]
	newfile2 = open(rootDir+'/'+'fileLocation', 'w')
	newfile2.write('\n'.join(files))
	newfile2.close()
	
def fastq_to_fasta(rootDir):
	print('Doing fastq to fasta conversion, wish me luck...')
	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0 and len(fileList) >= 1:
			folder = dirName.split('/')[-1]
			if folder.startswith('barcode') or folder.startswith('unclassified'):
				if len(list(Path(dirName).glob('*_renamed.fq'))) == 1 and len(list(Path(dirName).glob('*_renamed.fasta'))) == 0 :
					fqFile = str(list(Path(dirName).glob('*_renamed.fq'))[0])
					faFile = fqFile.replace('_renamed.fq','_renamed.fasta')
					#xxx = (r'">%s\n"')
					#equal = '%4 =='
					#cmd = "awk '{if(NR %s 1) {printf(%s,substr($0,2));} else if(NR %s 2) print;}' %s > %s"
					#os.system(cmd %(equal,xxx,equal,fqFile,faFile))
					cmd = "awk '{if(NR%4==1) {printf(\">%s\\n\",substr($0,2));} else if(NR%4==2) print;}'"
					os.system( cmd + ' %s > %s' %(fqFile,faFile) )
				elif len(list(Path(dirName).glob('*_renamed.fq'))) == 1 and len(list(Path(dirName).glob('*_renamed.fasta'))) == 1:
					continue
				else:
					print('Problematic folder = ' + dirName)
					
def rename_reads2(rootDir):
	dirnames = []
	seqPat = re.compile(r'sample\w*id=(\w+)')
	sampleid_dic = {'Cc_SE_1H':'C1H','Cc_SE_2H':'C2H','Cc_SE_3H':'C3H','Cc_SE_4H':'C4H','Cc_sE_5H':'C5H','Cc_sE_6H':'C6H'}
	for dirName, subdirList, fileList in os.walk(rootDir):
		dirnames += subdirList
		if len(subdirList) == 0 and len(fileList) >= 1:
			'''
			if len(list(Path(dirName).glob('*_renamed.fq'))) == 1 and len(list(Path(dirName).glob('*_renamed.fasta'))) == 1 and 'repeat_sequencing' not in dirnames:
				x = 0
				faFile = str(list(Path(dirName).glob('*_renamed.fasta'))[0])
				fqFile = str(list(Path(dirName).glob('*_renamed.fq'))[0])
				newFaFile = open(faFile.replace('_renamed.fasta','_renamed2.fasta'),'w')				
				length = int(os.popen('wc -l %s' %(faFile)).read().split()[0]) #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1

				file1=open(faFile,'r')
				print('Processing file ' + faFile)
				bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
				bar.start()
				while x < length: #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
					name=file1.readline()
					sequence=file1.readline()
					name_ed = name.split('_')[0]
					try:
						BAR = name.split()[0].split('_')[1]
					except:
						if name.split()[0].endswith('unclassified'):
							BAR = 'unclas'
							name_ed = name.split()[0].split('_')[0].replace('unclassified','')
					#BAR = name.split()[0].split('_')[1]
					sampleid = seqPat.search(name).group(1) #name.split()[2].replace('sampleid=','')
					if sampleid in sampleid_dic.keys():
						sampleid = sampleid_dic[sampleid]
					newFaFile.write(name_ed+'_'+sampleid+'_'+BAR+'\n'+sequence.strip()+'\n')		
					x+=2
					bar.update(x)
				bar.finish()
				x = 0
				length = int(os.popen('wc -l %s' %(fqFile)).read().split()[0])
				newFqFile = open(fqFile.replace('_renamed.fq','_renamed2.fq'),'w')
				file1=open(fqFile,'r')
				print('Processing file ' + fqFile)
				bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
				bar.start()
				while x < length: #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
					name=file1.readline()
					sequence=file1.readline()
					blah = file1.readline()
					score = file1.readline()
					name_ed = name.split('_')[0]
					try:
						BAR = name.split()[0].split('_')[1]
					except:
						if name.split()[0].endswith('unclassified'):
							BAR = 'unclas'
							name_ed = name.split()[0].split('_')[0].replace('unclassified','')
					#BAR = name.split()[0].split('_')[1]
					sampleid = seqPat.search(name).group(1) #name.split()[2].replace('sampleid=','')
					if sampleid in sampleid_dic.keys():
						sampleid = sampleid_dic[sampleid]
					newFqFile.write(name_ed+'_'+sampleid+'_'+BAR+'\n'+sequence.strip()+'\n'+blah.strip()+'\n'+score.strip()+'\n')		
					x+=4
					bar.update(x)
				bar.finish()
				x = 0
			'''
			if len(list(Path(dirName).glob('*_renamed.fq'))) == 1 and len(list(Path(dirName).glob('*_renamed.fasta'))) == 1: # and 'repeat_sequencing' in dirnames:
				x = 0
				print('yes')
				coding = {'CB1_BAR01':'C2H_BAR11_rpt',
					   'CB1_BAR02':'C2H_BAR12_rpt',
					   'CB1_BAR03':'C3H_BAR03_rpt',
					   'CB1_BAR04':'C3H_BAR11_rpt',
					   'CB1_BAR05':'C3H_BAR12_rpt',
					   'CB1_BAR06':'C5H_BAR10_rpt',
					   'CB1_BAR07':'C5H_BAR12_rpt',
					   'CB1_BAR08':'C6H_BAR01_rpt',
					   'CB1_BAR09':'C8H_BAR09_rpt','CB1_BAR10':'CB1_BAR10_rpt','CB1_BAR11':'CB1_BAR11_rpt','CB1_BAR12':'CB1_BAR12_rpt','CB1_unclas':'CB1_BAR09_unclas_rpt',
					   ##
						  'CB1r_BAR02':'C2H_BAR12_rpt',
					   'CB1r_BAR03':'C3H_BAR03_rpt',
					   'CB1r_BAR04':'C3H_BAR11_rpt',
					   'CB1r_BAR05':'C3H_BAR12_rpt',
					   'CB1r_BAR06':'C5H_BAR10_rpt',
					   'CB1r_BAR07':'C5H_BAR12_rpt',
					   'CB1r_BAR08':'C6H_BAR01_rpt',
					   'CB1r_BAR09':'C8H_BAR09_rpt','CB1r_BAR10':'CB1_BAR10_rpt','CB1r_BAR11':'CB1_BAR11_rpt','CB1r_BAR12':'CB1_BAR12_rpt','CB1r_unclas':'CB1_BAR09_unclas_rpt',
					   ##
					   'CB2r_BAR01':'C3H_BAR10_rpt',
					   'CB2r_BAR02':'C5H_BAR07_rpt',
					   'CB2r_BAR03':'C6H_BAR08_rpt',
					   'CB2r_BAR04':'C7H_BAR10_rpt',
					   'CB2r_BAR05':'C9H_BAR03_rpt',
					   'CB2r_BAR06':'C9H_BAR06_rpt',
					   'CB2r_BAR07':'C9H_BAR12_rpt',
					   'CB2r_BAR08':'C11H_BAR10_rpt',
					   'CB2r_BAR09':'C12H_BAR06_rpt', 'CB2r_BAR10':'CB2r_BAR10_rpt','CB2r_BAR11':'CB2r_BAR11_rpt','CB2r_BAR12':'CB2r_BAR12_rpt','CB2r_unclas':'CB2r_BAR09_unclas_rpt',
					   ##
						  'CB2_BAR01':'C3H_BAR10_rpt',
					   'CB2_BAR02':'C5H_BAR07_rpt',
					   'CB2_BAR03':'C6H_BAR08_rpt',
					   'CB2_BAR04':'C7H_BAR10_rpt',
					   'CB2_BAR05':'C9H_BAR03_rpt',
					   'CB2_BAR06':'C9H_BAR06_rpt',
					   'CB2_BAR07':'C9H_BAR12_rpt',
					   'CB2_BAR08':'C11H_BAR10_rpt',
					   'CB2_BAR09':'C12H_BAR06_rpt', 'CB2_BAR10':'CB2_BAR10_rpt','CB2_BAR11':'CB2_BAR11_rpt','CB2_BAR12':'CB2_BAR12_rpt','CB2_unclas':'CB2_BAR09_unclas_rpt',
						  #
					   'C15Hrpt2_BAR01':'C15H_BAR01_rpt',
					   'C15Hrpt2_BAR02':'C15H_BAR02_rpt',
					   'C15Hrpt2_BAR03':'C15H_BAR04_rpt',
					   'C15Hrpt2_BAR04':'C15H_BAR05_rpt',
					   'C15Hrpt2_BAR05':'C15H_BAR06_rpt',
					   'C15Hrpt2_BAR06':'C15H_BAR07_rpt',
					   'C15Hrpt2_BAR07':'C15H_BAR09_rpt',
					   'C15Hrpt2_BAR08':'C15H_BAR10_rpt',
					   'C15Hrpt2_BAR09':'C15H_BAR11_rpt',
					   'C15Hrpt2_BAR10':'C15H_BAR12_rpt',
					   'C15Hrpt2_BAR11':'C15H_BAR03_rpt',
					   'C15Hrpt2_BAR12':'C15H_BAR08_rpt',
					   'C15Hrpt2_unclas':'C15H_unclas_rpt',
					   ##
					   'C13Hr2_BAR01':'C13H_BAR01_rpt',
					   'C13Hr2_BAR02':'C13H_BAR02_rpt',
					   'C13Hr2_BAR03':'C13H_BAR03_rpt',
					   'C13Hr2_BAR04':'C13H_BAR04_rpt',
					   'C13Hr2_BAR05':'C13H_BAR05_rpt',
					   'C13Hr2_BAR06':'C13H_BAR06_rpt',
					   'C13Hr2_BAR07':'C13H_BAR07_rpt',
					   'C13Hr2_BAR08':'C13H_BAR08_rpt',
					   'C13Hr2_BAR09':'C13H_BAR09_rpt',
					   'C13Hr2_BAR10':'C13H_BAR10_rpt',
					   'C13Hr2_BAR11':'C13H_BAR11_rpt',
					   'C13Hr2_BAR12':'C13H_BAR12_rpt',
					   'C13Hr2_unclas':'C13H_unclas_rpt'}
				faFile = str(list(Path(dirName).glob('*_renamed.fasta'))[0])
				fqFile = str(list(Path(dirName).glob('*_renamed.fq'))[0])
				newFaFile = open(faFile.replace('_renamed.fasta','_renamed2.fasta'),'w')				
				length = int(os.popen('wc -l %s' %(faFile)).read().split()[0]) #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1

				file1=open(faFile,'r')
				print('Processing file ' + faFile)
				bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
				bar.start()
				while x < length: #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
					name=file1.readline()
					sequence=file1.readline()
					name_ed = name.split('_')[0]
					try:
						BAR = name.split()[0].split('_')[1]
					except:
						if name.split()[0].endswith('unclassified'):
							BAR = 'unclas'
							name_ed = name.split()[0].split('_')[0].replace('unclassified','')
					#BAR = name.split()[0].split('_')[1]
					sampleid = seqPat.search(name).group(1) #name.split()[2].replace('sampleid=','')
					if sampleid in sampleid_dic.keys():
						sampleid = sampleid_dic[sampleid]
					if sampleid+'_'+BAR in coding.keys():
						newsuffix = coding[sampleid+'_'+BAR]
						newFaFile.write(name_ed+'_'+newsuffix+'\n'+sequence.strip()+'\n')
					x+=2
					bar.update(x)
				bar.finish()
				x = 0
				length = int(os.popen('wc -l %s' %(fqFile)).read().split()[0])
				newFqFile = open(fqFile.replace('_renamed.fq','_renamed2.fq'),'w')
				file1=open(fqFile,'r')
				print('Processing file ' + fqFile)
				bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
				bar.start()
				while x < length: #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
					name=file1.readline()
					sequence=file1.readline()
					blah = file1.readline()
					score = file1.readline()
					name_ed = name.split('_')[0]
					try:
						BAR = name.split()[0].split('_')[1]
					except:
						if name.split()[0].endswith('unclassified'):
							BAR = 'unclas'
							name_ed = name.split()[0].split('_')[0].replace('unclassified','')
					#BAR = name.split()[0].split('_')[1]
					sampleid = seqPat.search(name).group(1) #name.split()[2].replace('sampleid=','')
					if sampleid in sampleid_dic.keys():
						sampleid = sampleid_dic[sampleid]
					if sampleid+'_'+BAR in coding.keys():
						newsuffix = coding[sampleid+'_'+BAR]
						newFqFile.write(name_ed+'_'+newsuffix+'\n'+sequence.strip()+'\n'+blah.strip()+'\n'+score.strip()+'\n')		
					x+=4
					bar.update(x)
				bar.finish()
				x = 0
			dirnames = []
			
def collect_fastq_files(rootDir):
	fqFilesWritten = []
	os.system("rm /home/banthony/scratch/reads/nanopore/capitata/all_cap_pass.fastq")
	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(subdirList) == 0 and len(fileList) >= 1 and 'repeat_sequencing' not in dirName.split('/'):			
			if len(list(Path(dirName).glob('*_renamed2.fq'))) == 1:
				fqFile = str(list(Path(dirName).glob('*_renamed2.fq'))[0])
				combine = "cat %s >> /home/banthony/scratch/reads/nanopore/capitata/all_cap_pass.fastq" %fqFile
				os.system(combine)
				fqFilesWritten += [fqFile]
	newFqFile = open('/home/banthony/scratch/reads/nanopore/capitata/fqFilesconcatenated','w')
	newFqFile.write('\n'.join(fqFilesWritten))
	newFqFile.close()

def create_links(rootDir):
	for dirName, subdirList, fileList in os.walk(rootDir):
		if str(dirName).split('/')[-1] == 'tofu2' and len(list(Path(dirName).glob('*min_fl_5.filtered.gff'))) >= 1:
			print("Working on %s" %dirName)
			grouptxt = str(list(Path(dirName).glob('*tofu.collapsed.group.txt'))[0]).split("/")[-1]
			gff = str(list(Path(dirName).glob('*tofu.collapsed.min_fl_5.filtered.gff'))[0]).split("/")[-1]
			counttxt = str(list(Path(dirName).glob('*tofu.collapsed.min_fl_5.filtered.abundance.txt'))[0]).split("/")[-1]
			repfq = str(list(Path(dirName).glob('*tofu.collapsed.min_fl_5.filtered.rep.fa'))[0]).split("/")[-1]
			#dirName = str(list(Path(dirName)[0]))
			#dirName = '/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/C1H/tofu3'
			os.system("cd %s && ln -s %s %s" %(str(dirName),grouptxt,'touse.group.txt'))
			os.system("cd %s && ln -s %s %s" %(str(dirName),gff,'touse.gff'))
			os.system("cd %s && ln -s %s %s" %(str(dirName),counttxt,'touse.count.txt'))
			os.system("cd %s && ln -s %s %s" %(str(dirName),repfq,'touse.rep.fq'))
	
def rename_reads3(fastafile):
	x = 0
	readnames = []
	length = int(os.popen('wc -l %s' %(fastafile)).read().split()[0])
	#newFaFile = open(fastafile.replace('correctedReads.fasta','correctedReads2.fasta'),'w')
	file1=open(fastafile,'r')
	print('Processing file ' + fastafile)
	bar = progressbar.ProgressBar(maxval=length+1, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	while x < length: #Because wc counts '\n' in the file sometimes if raises a ValueError: Value out of range. So u can add or remove 1
		name=file1.readline()
		sequence=file1.readline()
		#newFaFile.write(name.strip().split()[0]+'_corr\n'+sequence.strip()+'\n')
		readnames += [name.strip().split()[0].replace('_corr','')]		
		x+=2
		bar.update(x)
	bar.finish()
	#newFaFile.close()
	#os.system(" rm %s " %fastafile)
	return(readnames)

def collect_uncorrected_reads(rootDir):
	a=b=c=d=e = 0
	for dirName, subdirList, fileList in os.walk(rootDir):
		if len(list(Path(dirName).glob('*correctedReads2.fasta'))) == 1 and 'canu-logs' in subdirList: #len(list(Path(dirName).glob('*correctedReads.fasta.gz'))) == 1 and 'canu-logs' in subdirList:
			filename = list(Path(dirName).glob('*correctedReads2.fasta'))[0] #filename = list(Path(dirName).glob('*correctedReads.fasta.gz'))[0]
			print('Working on ' + str(filename))
			#os.system(" cd %s && gunzip -k %s " %(dirName,os.path.basename(filename)))
			print('File unzipped... Now renaming reads')
			readnames = rename_reads3(str(filename).replace('.gz',''))
			print('Reads renamed... Now getting original fasta dic')
			newDir = dirName.strip( dirName.split('/')[-1] ).replace('canu','pychoper')
			orgnFasta = list(Path(newDir).glob('*H_pass_py.all.fasta'))[0]
			
			orgnFasta_dic = find_fasta_seqs.fasta_to_dict(orgnFasta)
			#print(list(orgnFasta_dic.keys())[0])
			total_reads = len(orgnFasta_dic.keys())
			print('Querrying original fasta for corrected reads...')
			for read in readnames:
				if read.strip('>')+'\n' in orgnFasta_dic.keys():
					del orgnFasta_dic[read.strip('>')+'\n']
					a += 1
				else:
					b += 1
			orgnFastaNew = str(orgnFasta).replace('H_pass_py.all.fasta','_pass.uncorrected.fasta').replace('/pychoper/','/canu/')
			newFaFile = open(orgnFastaNew,'w')
			print('Querrying done, now writing new fasta with uncorrected reads.')
			for readx,sequencex in orgnFasta_dic.items():
				newFaFile.write('>'+readx+sequencex)
				c += 1
			newFaFile.close()
			print('Total reads = ' + str(total_reads) + '\nCor reads = ' + str(a) + '\nUn_Cor reads = ' + str(c) + '\nCor reads not found = '+str(b))
			print('Percentage Cor reads = ' + str(round(100*a/total_reads,2)))
			orgnFasta_dic=readnames = {}
			print('Now combining corrected and uncorrected reads...')
			combined_reads = str(orgnFasta).replace('H_pass_py.all.fasta','_pass.corNuncor_reads.fasta').replace('/pychoper/','/canu/')
			os.system(" cat %s %s > %s " %(str(str(filename).replace('.gz','')).replace('correctedReads.fasta','correctedReads2.fasta'),orgnFastaNew,combined_reads))
			print('Getting full length fasta file....')
			full_length = str(orgnFasta).replace('H_pass_py.all.fasta','_pass.FL_canu.fasta').replace('/pychoper/','/canu/')
			newFL = open(full_length, 'w')			
			for readY, sequenceY in find_fasta_seqs.fasta_to_dict(combined_reads).items():
				if '_FL' in readY:
					newFL.write(readY+sequenceY)
			newFL.close()
			print('Yeai, .... all done, moving on if theres another sample \n\n')
			a=b=c=d=e = 0

#count_files_and_sizes('/home/banthony/scratch')
#combine_fastq_files('/home/banthony/scratch/reads/nanopore/capitata/C011_9_8')
#combine_fastq_files('/home/banthony/scratch/reads/nanopore/capitata/C011_9_15/C1H/pass/')
#rename_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_15/C1H/fail/')
#rename_reads('/home/banthony/scratch/reads/nanopore/capitata/repeat_sequencing/')
#fastq_to_fasta('/home/banthony/scratch/reads/nanopore/capitata/repeat_sequencing/')
#rename_reads2('/home/banthony/scratch/reads/nanopore/capitata/Test/repeat_sequencing/barcode11')
#rename_reads2('/home/banthony/scratch/reads/nanopore/capitata/repeat_sequencing')
#rename_reads2('/home/banthony/scratch/reads/nanopore/capitata/repeat_sequencing/C011_09_B1')
#rename_reads2('/home/banthony/scratch/reads/nanopore/capitata/repeat_sequencing/C011_09_B2')
#collect_fastq_files('/home/banthony/scratch/reads/nanopore/capitata')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/Test')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_1')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_2')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_3')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_4')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_5')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_6')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_7')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_8')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_9')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_10')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_11')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_12')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_13')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_14')
#collect_uncorrected_reads('/home/banthony/scratch/reads/nanopore/capitata/C011_9_15')
create_links('/home/banthony/scratch/analysis/nanopore/capitata/rna_seq/single_embryo/')