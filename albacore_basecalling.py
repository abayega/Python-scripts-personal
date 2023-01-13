#!/usr/bin/env python3

import os, re
import progressbar

def control_jobs(source,flowcell,kit,batch,threads):
	filess = os.listdir(source)
	filess = ([x for x in filess if os.path.isdir(source+'/'+x)])
	
	for f in filess:
		if f != 'pass' and f != 'fail':
			destination = source+'/'+f+'/basecalled'
			destination1 = source+'/'+f+'/jobs'
			sourcex = source+'/'+f+'/fast5'
			#print(destination1)
			if not os.path.isdir(destination):
				os.makedirs(destination)
			if not os.path.isdir(destination1):
				os.makedirs(destination1)
			files = os.listdir(sourcex)
			nfiles = len(files)
			flowcell = flowcell #'FLO-MIN106'
			kit = kit #'SQK-LSK108'
			threads = str(threads)
			brk = int(batch)
			python3_path = "read_fast5_basecaller.py" # or ~/.local/bin/read_fast5_basecaller.py
			n=0
			x=0
			step = nfiles/brk
			if '.' in str(step):
				k = int(str(step).split('.')[0])+1
			else:
				k = step
			for i in range(0,k):
				newfile = open(destination1+'/albacore_%s' %(i)+'.py', 'w')
				if len(files[n:]) >= brk:
					job_def = '#!/usr/bin/env python3\n\n'+'import os, re, progressbar\n\n'+'def submit_jobs3(source,destination):\n'+'\tfiles = sorted(os.listdir(source))\n'+'\tflowcell = %s\n' %('"'+flowcell+'"')+'\tkit = %s\n' %('"'+kit+'"')+'\tthreads = %s\n' %('"'+threads+'"')+'\tpython3_path = %s\n' %('"'+python3_path+'"')+'\tn=0\n'+'\tprint("files[%s:%s]")\n' %(n,n+brk)+\
					'\tprint(files[%s:%s])\n' %(n,n+brk)+"\t#bar = progressbar.ProgressBar(maxval=len(files[%s:%s]), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])\n" %(n,n+brk)+'\t#bar.start()\n'+\
					'\tfor i in files[%s:%s]:\n' %(n,n+brk) +'\t\tinput_dir = source+"/"+i\n'+'\t\tsave_path = destination+"/"+i\n'+'\t\tif not os.path.exists(save_path):\n'+'\t\t\tos.makedirs(save_path)\n'+\
					'\t\tprint("File " + str(n+1) + " of " + str(%s))\n' %(len(files[n:n+brk]))+'\t\tos.system("%s -r --flowcell %s --kit %s --input %s --save_path %s --worker_threads %s -o fastq" %(python3_path,flowcell,kit,input_dir,save_path,int(threads)))\n'+'\t\tn+=1\n'+'\t\t#bar.update(n)\n'+'\t#bar.finish()\n\n'+'submit_jobs3("%s","%s")' %(sourcex,destination)
				else:
					job_def = '#!/usr/bin/env python3\n\n'+'import os, re, progressbar\n\n'+'def submit_jobs3(source,destination):\n'+'\tfiles = sorted(os.listdir(source))\n'+'\tflowcell = %s\n' %('"'+flowcell+'"')+'\tkit = %s\n' %('"'+kit+'"')+'\tthreads = %s\n' %('"'+threads+'"')+'\tpython3_path = %s\n' %('"'+python3_path+'"')+'\tn=0\n'+'\tprint("files[%s:%s]")\n' %(n,nfiles)+\
					'\tprint(files[%s:%s])\n' %(n,n+brk)+"\t#bar = progressbar.ProgressBar(maxval=len(files[%s:%s]), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])\n" %(n,nfiles)+'\t#bar.start()\n'+\
					'\tfor i in files[%s:]:\n' %(n) +'\t\tinput_dir = source+"/"+i\n'+'\t\tsave_path = destination+"/"+i\n'+'\t\tif not os.path.exists(save_path):\n'+ '\t\t\tos.makedirs(save_path)\n'+\
					'\t\tprint("File " + str(n+1) + " of " + str(%s))\n' %(len(files[n:]))+'\t\tos.system("%s -r --flowcell %s --kit %s --input %s --save_path %s --worker_threads %s -o fastq" %(python3_path,flowcell,kit,input_dir,save_path,int(threads)))\n'+'\t\tn+=1\n'+'\t\t#bar.update(n)\n'+'\t#bar.finish()\n\n'+'submit_jobs3("%s","%s")' %(sourcex,destination)			
				newfile.write(job_def)
				newfile.close()			
				n+=brk

			#now submit these jobs
			files = ([x for x in os.listdir(destination1) if x.endswith('.py')])
			for i in files:
				newfile = open(destination1+'/'+ i.replace('py','sh'), 'w')
				job_def = '#!/bin/bash\n'+'#PBS -A vsy-012-01\n'+'#PBS -N %s\n' %(i.replace('.py',''))+'#PBS -l walltime=24:00:00\n'+ '#PBS -l nodes=1:ppn=1\n'+'#PBS -q qwork\n'+'#PBS -r n\n\n'+'module load python64/3.5.2\n\n'+\
				'export PATH=~/.local/bin:$PATH #This is the PATH where pip3 stores site executables\n\n'+'dirc=%s\n\n' %(destination1)+'cd $dirc\n\n'+'echo "starting albacore"\n'+'python3 %s' %(destination1+'/'+ i) + '\necho "albacore finished"'
				newfile.write(job_def)
				newfile.close()
				#os.system('qsub %s' %(destination1+'/'+ i.replace('py','sh')))		

def control_jobsx(source,destination,destination1):
	files = os.listdir(source)
	nfiles = len(files)
	brk = 90
	n=0
	x=0
	step = nfiles/brk
	flowcell = 'FLO-MIN106'
	kit = 'SQK-LSK108'
	if '.' in str(step):
		k = int(str(step).split('.')[0])+1
	else:
		k = step
	for i in range(0,k):
		newfile = open(destination1+'/albacore_%s' %(i)+'.py', 'w')
		if len(files[n:]) >= brk:
			job_def = '#!/usr/bin/env python3\n\n'+'import os, re, progressbar\n\n'+'def submit_jobs3(source,destination):\n'+'\tfiles = sorted(os.listdir(source))\n'+'\tflowcell = %s\n' %('"'+flowcell+'"')+'\tkit = %s\n' %('"'+kit+'"')+'\tn=0\n'+'\tprint("files[%s:%s]")\n' %(n,n+brk)+\
			'\tprint(files[%s:%s])\n' %(n,n+brk)+"\t#bar = progressbar.ProgressBar(maxval=len(files[%s:%s]), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])\n" %(n,n+brk)+'\t#bar.start()\n'+\
			'\tfor i in files[%s:%s]:\n' %(n,n+brk) +'\t\tinput_dir = source+"/"+i\n'+'\t\tsave_path = destination+"/"+i\n'+'\t\tif not os.path.exists(save_path):\n'+'\t\t\tos.makedirs(save_path)\n'+\
			'\t\tprint("File " + str(n+1) + " of " + str(%s))\n' %(len(files[n:n+brk]))+'\t\tos.system("~/.local/bin/read_fast5_basecaller.py -r --flowcell %s --kit %s --input %s --save_path %s --worker_threads 23 -o fastq" %(flowcell,kit,input_dir,save_path))\n'+'\t\tn+=1\n'+'\t\t#bar.update(n)\n'+'\t#bar.finish()\n\n'+'submit_jobs3("%s","%s")' %(source,destination)
		else:
			job_def = '#!/usr/bin/env python3\n\n'+'import os, re, progressbar\n\n'+'def submit_jobs3(source,destination):\n'+'\tfiles = sorted(os.listdir(source))\n'+'\tflowcell = %s\n' %('"'+flowcell+'"')+'\tkit = %s\n' %('"'+kit+'"')+'\tn=0\n'+'\tprint("files[%s:%s]")\n' %(n,nfiles)+\
			'\tprint(files[%s:%s])\n' %(n,n+brk)+"\t#bar = progressbar.ProgressBar(maxval=len(files[%s:%s]), widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])\n" %(n,nfiles)+'\t#bar.start()\n'+\
			'\tfor i in files[%s:]:\n' %(n) +'\t\tinput_dir = source+"/"+i\n'+'\t\tsave_path = destination+"/"+i\n'+'\t\tif not os.path.exists(save_path):\n'+ '\t\t\tos.makedirs(save_path)\n'+\
			'\t\tprint("File " + str(n+1) + " of " + str(%s))\n' %(len(files[n:]))+'\t\tos.system("~/.local/bin/read_fast5_basecaller.py -r --flowcell %s --kit %s --input %s --save_path %s --worker_threads 23 -o fastq" %(flowcell,kit,input_dir,save_path))\n'+'\t\tn+=1\n'+'\t\t#bar.update(n)\n'+'\t#bar.finish()\n\n'+'submit_jobs3("%s","%s")' %(source,destination)			
		newfile.write(job_def)
		newfile.close()			
		n+=brk

	#now submit these jobs
	files = ([x for x in os.listdir(destination1) if x.endswith('.py')])
	for i in files:
		newfile = open(destination1+'/'+ i.replace('py','sh'), 'w')
		job_def = '#!/bin/bash\n'+'#PBS -A vsy-012-01\n'+'#PBS -N %s\n' %(i.replace('.py',''))+'#PBS -l walltime=24:00:00\n'+ '#PBS -l nodes=1:ppn=1\n'+'#PBS -q qwork\n'+'#PBS -r n\n\n'+'module load python64/3.5.2\n\n'+\
		'export PATH=~/.local/bin:$PATH #This is the PATH where pip3 stores site executables\n\n'+'dirc=%s\n\n' %(destination1)+'cd $dirc\n\n'+'echo "starting albacore"\n'+'python3 %s' %(destination1+'/'+ i) + '\necho "albacore finished"'
		newfile.write(job_def)
		newfile.close()
		#os.system('qsub %s' %(destination1+'/'+ i.replace('py','sh')))		
		
def collect_fastq(source,destination_file,pass_or_fail):
	files = os.listdir(source)
	files = ([x for x in files if os.path.isdir(source+'/'+x)])
	n=0
	m=len(files)
	print('Source directory has ' + str(m) + ' directories')
	k = 0
	bar = progressbar.ProgressBar(maxval=m, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for file in files:
		if 'basecalled' in os.listdir(source+'/'+file):
			basecalled_dir = source+'/'+file+'/'+'basecalled/'
			files2 = os.listdir(basecalled_dir)
			for x in files2:
				if x.isdecimal():
					path = basecalled_dir+x+'/workspace/'+pass_or_fail
					os.system("cat %s/* >> %s" %(path,destination_file))
					k += 1
		n+=1
		bar.update(n)
	bar.finish()
	print('A total of ' + str(k) + ' files were concatenated')
	
def collect_fastq_barcoded(source,destiation_file_unclassified,destiation_file_barcode02,destiation_file_barcode08,pass_or_fail):
	files = os.listdir(source)
	files = ([x for x in files if os.path.isdir(source+'/'+x)])
	n=0
	m=len(files)
	print('Source directory has ' + str(m) + ' directories')
	k = 0
	bar = progressbar.ProgressBar(maxval=m, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for file in files:
		if 'basecalled' in os.listdir(source+'/'+file):
			basecalled_dir = source+'/'+file+'/'+'basecalled/'
			files2 = os.listdir(basecalled_dir)
			for x in files2:
				if x.isdecimal():					
					path = basecalled_dir+x+'/workspace/'+pass_or_fail
					dir_path = os.listdir(path)
					if 'unclassified' in dir_path:						
						newpath = basecalled_dir+x+'/workspace/'+pass_or_fail+'/unclassified'
						os.system("cat %s/* >> %s" %(newpath,destination_file_unclassified))
					if 'barcode02' in dir_path:
						newpath = basecalled_dir+x+'/workspace/'+pass_or_fail+'/barcode02'
						os.system("cat %s/* >> %s" %(newpath,destination_file_barcode02))
					if 'barcode08' in dir_path:
						newpath = basecalled_dir+x+'/workspace/'+pass_or_fail+'/barcode08'
						os.system("cat %s/* >> %s" %(newpath,destination_file_barcode08))					
					k += 1
		n+=1
		bar.update(n)
	bar.finish()
	print('A total of ' + str(k) + ' files were concatenated')
	
def collect_seqstats(source,destination_file):
	files = os.listdir(source)
	files = ([x for x in files if os.path.isdir(source+'/'+x)])
	n=0
	m=len(files)
	print('Source directory has ' + str(m) + ' directories')
	k = 0
	bar = progressbar.ProgressBar(maxval=m, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for file in files:
		if 'basecalled' in os.listdir(source+'/'+file):
			basecalled_dir = source+'/'+file+'/'+'basecalled/'
			files2 = os.listdir(basecalled_dir)
			for x in files2:
				if x.isdecimal():
					path = basecalled_dir+x+'/sequencing_summary.txt'
					os.system("cat %s >> %s" %(path,destination_file))
					k += 1
		n+=1
		bar.update(n)
	bar.finish()
	print('A total of ' + str(k) + ' files were concatenated')
	
def remove_basecalled_fast5(source):
	files = os.listdir(source)
	files = ([x for x in files if os.path.isdir(source+'/'+x)])
	n=0
	m=len(files)
	print('Source directory has ' + str(m) + ' directories')
	k = 0
	k1 = 0
	bar = progressbar.ProgressBar(maxval=m, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
	bar.start()
	for file in files:
		if 'basecalled' in os.listdir(source+'/'+file):
			basecalled_dir = source+'/'+file+'/'+'basecalled/'
			os.system("rm -r %s" %(basecalled_dir))
			k += 1
		if 'fast5' in os.listdir(source+'/'+file):
			fast5_dir = source+'/'+file+'/'+'fast5/'
			os.system("rm -r %s" %(fast5_dir))
			k1 += 1
		n+=1
		bar.update(n)
	bar.finish()
	print('A total of ' + str(k) + ' basecalled directories and ' + str(k1) + ' fast5 directories were deleted')
			
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/reads/nanopore/river_sequencing/C008_07_290118'
dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2017/ioannisr/banthony/course_2018/s4'
#dr = '/mnt/parallel_scratch_mp2_wipe_on_december_2017/ioannisr/banthony/new_course_files_5/new_file_2/s12'

destination_file = dr+'/pass/s4_basecalled_pass.fastq'
#destination_file_unclassified = dr+'/fail/unclassified/ONT_River_C008_05_fail_unclassfied.fastq'
#destination_file_barcode02 = dr+'/fail/barcode02/ONT_River_C008_05_fail_barcode02.fastq'
#destination_file_barcode08 = dr+'/fail/barcode08/ONT_River_C008_05_fail_barcode08.fastq'


#control_jobs(dr,'FLO-MIN106','SQK-LSK108',10,23)
collect_fastq(dr,destination_file,'pass')
destination_file = dr+'/s4_sequencing_summary.txt'
collect_seqstats(dr,destination_file)
#remove_basecalled_fast5(dr)

#control_jobsx(dr+'/fast5',dr+'/basecalled',dr+'/jobs')
#collect_fastq_barcoded(dr,destination_file_unclassified,destination_file_barcode02,destination_file_barcode08,'fail')
#print(1)
#remove_basecalled_fast5(dr)
#print(2)
#dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/reads/nanopore/river_sequencing/C008_05_091117'
#remove_basecalled_fast5(dr)
#print(3)
#dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/reads/nanopore/olivefly/rna_seq/Bo_E_6H/C010_07_4_041117'
#remove_basecalled_fast5(dr)
#print(4)
#dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/reads/nanopore/olivefly/rna_seq/Bo_E_4H/C010_06_191017'
#remove_basecalled_fast5(dr)
#print(5)
#dr = '/mnt/parallel_scratch_mp2_wipe_on_august_2017/ioannisr/banthony/reads/nanopore/olivefly/rna_seq/Bo_E_2H/C010_09_141117'
#remove_basecalled_fast5(dr)
