#!/usr/bin/env python

import h5py, os, sys, stat, math

global newfile, sumbad, sumeve, sumreadln

def Get_CallStatsFromFile(fast5_path):
    global newfile, sumbad, sumeve, sumreadln
    sumbad = 0
    sumeve = 0
    sumreadln = 0
    'Parse the fast5 analysed file and return file name and fasta length.'
    fast5_filename = os.path.basename(fast5_path)
    run_number = fast5_filename.split('_')[-4]
    #run_number = '_'.join(fast5_filename.split('_')[-5:-3]) BA
    file_number = fast5_filename.split('_')[-2].replace('read', '')
    if os.stat(fast5_path)[stat.ST_SIZE] == 0:
        sys.stderr.write('Erro: analysed fast5 file has length of zero - ignoring ({0})\n'.format(fast5_path))
        return None
    try:
        hdf = h5py.File(fast5_path, 'r')
    except:
        #sys.stderr.write('Erro: Failed to open file with h5py.File for unknown reason - file must be corrupt ({0})\n'.format(fast5_path))
	sumbad = sumbad + 1
        return None
    try:
        read_number = [int(x.split('_')[1]) for x in hdf['Analyses/EventDetection_000/Reads'].keys() if x.startswith('Read_')][0]
    except:
        #sys.stderr.write('Erro: Failed to get read_number from EventDetection_000 section - file must be corrupt ({0})\n'.format(fast5_path))
   	sumeve = sumeve + 1
        return None

    try:
        # New data
        key = 'Analyses/Basecall_1D_000/Summary/basecall_1d_template'
	key2 = 'Analyses/EventDetection_000/Reads/Read_{0}'.format(read_number)
        #channel_number = int(hdf[key].attrs['channel_number'])
        #digitisation = float(hdf[key].attrs['digitisation'])
        #offset = float(hdf[key].attrs['offset'])
        #range = float(hdf[key].attrs['range'])
        #sampling_rate = float(hdf[key].attrs['sampling_rate'])
    except:
        # Old data
        key = 'UniqueGlobalKey/read_id'
        channel_number = int(hdf[key].attrs['channel_number'])
        read_number2 = int(hdf[key].attrs['read_number'])
        digitisation = 0.0
        offset = 0.0
        range = 0.0
        sampling_rate = _RS[fast5_filename] if fast5_filename in _RS.keys() else _RS[_RS.keys()[0]]

    try:
        readlength = hdf[key].attrs['sequence_length']
	read_id = hdf[key2].attrs['read_id']
    except:
        #sys.stderr.write('Erro: Failed to open file with h5py.File for unknown reason - file must be corrupt ({0})\n'.format(fast5_path))
        sumreadln = sumreadln + 1
        return None
    writ = fast5_filename + '\t' + read_id + '\t' + str(readlength) + '\n'
    
    newfile.write(str(writ))

def parsefile(directory):
    global newfile, sumbad, sumeve, sumreadln 
    files = os.listdir(directory)
    newfile = open('newfileFail.txt', 'w')    
    writ = 'fast5_filename' + '\t' + 'read_id' + '\t' + 'readlength' + '\n'
    
    newfile.write(str(writ))    
   
    files = os.listdir(directory)
    for afile in files:
	Get_CallStatsFromFile(os.path.join(directory,afile))

    newfile.close()
    print('Sum of unopenable files is ' + str(sumbad) + '\n' + 'sum of no event is ' + str(sumeve) + '\n' + 'Sum of no readlength is ' + str(sumreadln))



parsefile('/gs/project/wst-164-ab/anthony/Nanopore_data/olive_fly/C004_07B_4_310117/fail')
    
