#!/bin/env python3
### This script removes ISPCR in 2D reads after they are already demultiplexed.
### It assumes your 2D reads are in a file called 2D.fastq and will return two files called 2D_trimmed_l.fastq 2D_trimmed_l.fasta
### Read names will be appended with 'p's if adapters are found and removed on either end, 
### i.e. _n_p would mean that an adapter was found on the right but not the left end of the read
### You'll need to install the "editdistance" module for python

###This script is for fastq file processing

import editdistance
import sys, os
import progressbar

#path=sys.argv[1]
#input1='2D.fastq'
#output1='2D_trimmed_l'

'''Molecule architecture
5' AAGCAGTGGTATCAACGCAGAGTGGATTCTATCACGCGGG   AAAAAAA      GTTGCGTTGCATACTCTGCGTTGATACCACTGCTT 3' Forward (+) strand as sequenced by ONT
5' AAGCAGTGGTATCAACGCAGAGTATGCAACGCAAC        TTTTTTT CCCGCGTGATAGAATCCACTCTGCGTTGATACCACTGCTT 3' Reverse (-) strand as sequenced by ONT
'''

TSO='AAGCAGTGGTATCAACGCAGAGT' #23 forward molecule AAGCAGTGGTATCAACGCAGAGT GGATTCTATC ACGCGG xxxxxx AAAAAAAAAAAAAAA GAAGCGTTGCA TACTCTGCGTTGATACCACTGCTT
formk = 'GGATTCTATCACGCGG' #16
revmk = 'TTTTTTTTTTTTTTT' #15 TGCAACGCAAC' #15 TTTTTT' 
#revmk = 'TGCAACGCAACTTTTT' #15          
polyT = 'ACTCTGCGTTGATACCACTGCTT' #23 rev molecule AAGCAGTGGTATCAACGCAGAGT TGCAACGCAAC TTTTTTTTTTTTTTT xxxxxxx CCGCGT GATAGAATCC ACTCTGCGTTGATACCACTGCTT
formk2 = 'AAAAAAAAAAAAAAA' #15 # 'GAAGCGTTGCA' #10
#formk2 = 'AAAAAGAAGCGTTGCA'
revmk2 = 'CCGCGTGATAGAATCC' #16

Counter_dict={}
Counter_dict['GGATTCTATC']=0
Counter_dict['TGCAACGCAA']=0

#if not os.path.exists(path):
#  os.makedirs(path)

def reverse_complement(sequence):
  Seq=''  
  complement = {'a':'T','c':'G','g':'C','t':'A','n':'N','A':'T','C':'G','G':'C','T':'A','N':'N','-':'-'}
  for item in sequence[::-1]:
    Seq=Seq+complement[item]
  return Seq

def find_ispcr(infile):
  path = os.path.dirname(infile)
  outputa=os.path.basename(infile).replace('fastq','')+'trimmed_stranded'
  outputb=os.path.basename(infile).replace('fastq','')+'no_adapters'
  outputc=os.path.basename(infile).replace('fastq','')+'single_adapter'
  outputd=os.path.basename(infile).replace('fastq','')+'unknown'
  outq=open(path+'/'+outputa+'.fastq','w')
  #outa=open(path+'/'+outputa+'.fasta','w')
  #outa2=open(path+'/'+outputa+'.2fasta','w')
  outb=open(path+'/'+outputb+'.fastq','w')
  outc=open(path+'/'+outputc+'.fastq','w')
  outd=open(path+'/'+outputd+'.fastq','w')
  
  length=0
  for line in open(infile):
    length+=1

  #print(length)

  x=0

  file1=open(infile,'r')

  total=0
  forward=0
  reverse=0
  total_double=0
  difference=0
  neg_strand=0
  match_left=''
  match_right=''
  strand = 'un'
  rbackseen = 'nr'
  fbackseen = 'nf'
  wu1 = 0
  wu2 = 0
  wu3 = 0
  trim = 10  #use 10 to match the original Mandalorion algorithm or 20 to trim as much as possible
  bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
  bar.start()
  while x<length:
    match_left=''
    match_right=''
    strand = ''
    rbackseen = ''
    fbackseen = ''
    a=file1.readline().strip().split()[0]
    b=file1.readline().strip()
    c=file1.readline()
    d=file1.readline().strip()
    reverse_b=b[::-1]    
    length_seq=len(b)
    forward_trimm='-'
    reverse_trimm='-'
    total+=1
    mindist1=11
    mindist1f=5
    mindist1r=5
    for xx in range(0,100,1):
        forward_seq=b[xx:xx+23]        
        dist=editdistance.eval(forward_seq,TSO)
        #print(str(xx)+' '+forward_seq+' '+str(dist)+' '+TSO)
        if dist<mindist1:
            #print('\n'+str(xx)+' '+forward_seq+' '+str(dist)+' '+TSO)                                   
            for y in range((xx+23)-2,xx+43,1):
                distf = editdistance.eval(b[y:y+16],formk)
                distr = editdistance.eval(b[y:y+15],revmk)
                #print(str(y)+' '+b[y:y+16]+' '+str(editdistance.eval(b[y:y+16],formk))+' '+formk)
                if distf < mindist1f: #or b[y:y+4]=='TTTT':
                    strand = '+f'                
                    #match_left=b[y:y+16]
                    #match_left=b[y:y+4]
                    mindist1=dist
                    mindistf=distf
                    #forward_trimm=y+16+3
                    forward_trimm=y+trim
                    break
                elif distr < mindist1r:
                    strand = '-r'
                    #match_left=b[y:y+4]
                    mindist1=dist
                    mindist1r=distr
                    #forward_trimm=y+15+3
                    forward_trimm=y+trim
                    break
                #else:
                #    strand = 'un'
                #    match_left=''

    mindist2=11
    mindist2f=5
    mindist2r=5
    for xx in range(0,100,1):
        #print('moving to reverse')
        reverse_seq=reverse_b[xx:xx+23]        
        dist=editdistance.eval(reverse_seq,polyT[::-1])
        #print(str(xx)+' '+reverse_seq+' '+str(dist)+' '+polyT[::-1])
        if dist<mindist2:
            #print('found reverse match')
            #print('\n'+str(xx)+' '+reverse_seq+' '+str(dist)+' '+polyT[::-1])
            for y in range((xx+23)-2,xx+43,1):
                distf2 = editdistance.eval(reverse_b[y:y+15],formk2[::-1])
                distr2 = editdistance.eval(reverse_b[y:y+16],revmk2[::-1])
                #print(str(y)+' '+reverse_b[y:y+16]+' '+str(editdistance.eval(reverse_b[y:y+16],revmk2[::-1]))+' '+revmk2[::-1])
                if distf2 < mindist2f: # or reverse_b[y:y+4]=='TTTT':
                    fbackseen = 'fb'
                    #strand = '+f'
                    #Counter_dict[reverse_b[y:y+10]]+=1
                    #match_right=reverse_b[y:y+4]
                    mindist2=dist
                    mindist2f=distf2
                    #reverse_trimm=length_seq-(y+15+3)
                    reverse_trimm=length_seq-(y+trim)
                    break
                elif distr2 < mindist2r: # or reverse_b[y:y+4]=='TTTT':
                    #print('yes'+' '+str(distr2)+' '+str(mindist2r))
                    rbackseen = 'rb'
                    #strand = '-r'
                    #Counter_dict[reverse_b[y:y+10]]+=1
                    #match_right=reverse_b[y:y+4]
                    mindist2=dist
                    mindist2r=distr2
                    #reverse_trimm=length_seq-(y+16+3)
                    reverse_trimm=length_seq-(y+trim)
                    break
                #else:
                #    rbackseen = 'nr'
                #    fbackseen = 'nf'
                #    match_right=''

    x+=4
    
    add_left='n'
    add_right='n'
    if forward_trimm!='-':
        forward+=1
        add_left='p' 
    else:
        forward_trimm=0
    if reverse_trimm!='-':
        reverse+=1
        add_right='p'
    else:
        reverse_trimm=length_seq
    if add_right!='n' and add_left!='n':
        total_double+=1
        #if match_left!=match_right:
        if strand == '+f':
            difference+=1
        elif strand == '-r':
            neg_strand+=1
    #print(match_left+' '+match_right)
    
    if add_right == 'n' and add_left == 'p':
        wu1 += 1
    elif add_right == 'p' and add_left == 'n':
        wu2 += 1
    elif add_right == 'n' and add_left == 'n':
        wu3 += 1
    if add_right == 'p' and add_left == 'p' and strand == '+f' and fbackseen == 'fb':        
        outq.write(a+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+b+'\n'+c+d+'\n')
        #outa.write('>'+a[1:]+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+b+'\n')
    elif add_right == 'p' and add_left == 'p' and strand == '-r' and rbackseen == 'rb':
        reversed_seq = reverse_complement(b)
        reversed_q = d[::-1]
        outq.write(a+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+reversed_seq+'\n'+c+reversed_q+'\n')
        #outa.write('>'+a.strip()[1:]+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+reversed_seq+'\n')
    elif add_right != 'p' and add_left != 'p':
        #outb.write('>'+a.strip()[1:].split()[0]+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+b[forward_trimm:reverse_trimm]+'\n')
        outb.write(a+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+b+'\n'+c+d+'\n')
    elif add_right == 'p' or add_left == 'p' and strand == '+f' or fbackseen == 'fb':
        outc.write(a+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+b+'\n'+c+d+'\n')
    elif add_right == 'p' or add_left == 'p' and strand == '-r' or rbackseen == 'rb':
        reversed_seq = reverse_complement(b)
        reversed_q = d[::-1]
        outc.write(a+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+reversed_seq+'\n'+c+reversed_q+'\n')
    else:
        outd.write(a+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+b+'\n'+c+d+'\n')

    #elif add_right != 'p' and add_left != 'p':
    #    outb.write('>'+a.strip()[1:]+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+b[forward_trimm:reverse_trimm]+'\n')
    #else:
    #    outc.write('>'+a.strip()[1:]+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+b[forward_trimm:reverse_trimm]+'\n')

    bar.update(x)

    #outq.write(a.strip()+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+b[forward_trimm:reverse_trimm]+'\n'+c+d[forward_trimm:reverse_trimm]+'\n')
    #outa2.write('>'+a.strip()[1:]+'.'+strand+'.'+fbackseen+'.'+rbackseen+'.'+'_'+add_left+'_'+add_right+'\n'+b[forward_trimm:reverse_trimm]+'\n')
  #print('Total reads = '+ str(total) + '\nTotal full length = ' + str(total_double) + '\n%ge of full length reads = ' + str("{0:.2f}".format(round(((total_double/total)*100),2))) + '\nIgnore = ' + str("{0:.2f}".format(round(((difference/total_double)*100),2))) + '\nTotal reads with forward adapter = ' + str(forward) + '\nTotal reads with back adapter = ' + str(reverse) +  '\n%ge of reads with fwd adapter = ' + str("{0:.2f}".format(round(((forward/total)*100),2))) + '\n%ge of reads with back adapter = ' + str("{0:.2f}".format(round(((reverse/total)*100),2))) + '\nignore = ' + str(difference) + '\nignore = ' + str("{0:.2f}".format(round(((difference/total_double)*100),2))))
  bar.finish()
  u=0
  try:
    print('Total reads = ' + str(total))
  except:
    u += 1
  try:
    print('Total full length = ' + str(total_double))
  except:
    u += 1
  try:
    print('%ge of full length reads = ' + str("{0:.2f}".format(round(((total_double/total)*100),2)))) 
  except:
    u += 1  
  try:
    print('Total reads with forward adapter = ' + str(forward))
  except:
    u += 1
  try:
    print('Total reads with back adapter = ' + str(reverse))
  except:
    u += 1
  try:
    print('%ge of reads with fwd adapter = ' + str("{0:.2f}".format(round(((forward/total)*100),2))))
  except:
    u += 1
  try:
    print('%ge of reads with back adapter = ' + str("{0:.2f}".format(round(((reverse/total)*100),2))))
  except:
    u +=1
  try:
    print('Total +Ve strand reads = ' + str(difference))
  except:
    u +=1
  try:
    print('%ge of +Ve strand reads to reads with both adapters = ' + str("{0:.2f}".format(round(((difference/total_double)*100),2))))
  except:
    u += 1
  try:
    print('%ge of -Ve strand  reads to reads with both adapters = ' + str("{0:.2f}".format(round(((neg_strand/total_double)*100),2))))
  except:
    u += 1
  try:
    print('%ge of reads with only one adapter = ' + str("{0:.2f}".format(round((((wu1+wu2)/total)*100),2))))
  except:
    u +=1  
  try:
    #print('%ge of reads with only reverse adapter = ' + str("{0:.2f}".format(round(((wu2/total)*100),2))))
    print('%ge of reads without any adapter = ' + str("{0:.2f}".format(round(((wu3/total)*100),2))))
    #print(Counter_dict)  
  except:
    u +=1
  try:
    print('Number of bases trimmed from each end = ' + str(trim))
  except:
    u +=1
  print('Failed to print = '+str(u))

find_ispcr('/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/reads/nanopore/olivefly/rna_seq/Bo_all/pass/mandalorion/fastq/Bo_E_2H_C010_09_pass_edited.fastq')
#find_ispcr('/home/abayega/Desktop/Untitled_Folder/capitatata.fq')
#find_ispcr('/home/abayega/Desktop/Untitled_Folder/sequences2.txt')


















