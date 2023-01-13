#!/bin/bash
#PBS -l nodes=1:ppn=12,pmem=1700M,walltime=06:00:00
#PBS -N nameofscript
#PBS -V
#PBS -A wst-164-ab
#PBS -m ae 
#PBS -W umask=0002 
#PBS -j oe 
#PBS -q metaq 

#qsub -l walltime=01:00:00 -l nodes=1:ppn=4,pmem=2700m -N ${fn}_bamstats -q sw -A wst-164-ab -e $errors/${fn}_err -o $out_log/${fn}_out

cd ~/../../gs/project/wst-164-ab/spyros/polyAplusandminus/normoxia/polyAminusN1

module load OpenBLAS+LAPACK/0.2.12-openmp-gcc
module load raglab/zlib/1.2.8
#module load gcc/4.8.2
#module load gcc/4.9.1
module load GCC/5.3.0-2.26 #Important for lastdb used nanook
module load raglab/redis/2.8.18
module load raglab/trimmomatic/0.33
module load raglab/bwa/0.7.15
module load raglab/bowtie2/2.2.3


#q1=`ls *_1.fastq.gz|cut -d " " -f 1`
#q2=`ls *_2.fastq.gz|cut -d " " -f 1`

gzip -d WTCHG_52392_279_1.fastq.gz
gzip -d WTCHG_52392_279_2.fastq.gz

#a1=`ls *_1.fastq|cut -d " " -f 1`
#a2=`ls *_2.fastq|cut -d " " -f 1`

tophat --no-discordant --no-coverage-search --transcriptome-max-hits 20 --transcriptome-index transcriptome_folder/known --no-mixed -r 500 --mate-std-dev 500  --library-type fr-firststrand -o alignment -p 12 -G ../../commongtf/selected_Homo_sapiens.GRCh37.75.gtf ../../commongenome/hg19 WTCHG_52392_279_1.fastq WTCHG_52392_279_2.fastq

ssh bicyclez@guillimin.hpc.mcgill.ca
qsub -A wst-164-ab guillimin_submit_npscarf.sh
showq -u bicyclez #this will show what is going on with your submitted jobs
qsub -A wst-164-ab <name of the file with commands>
qstat -u bicyclez #for queued jobs
qdel <jobID> #to cancel job
exit #to exit interactive node
qsub -I -A wst-164-ab -d `pwd` -l walltime=16:00:00 -q metaq -l nodes=1:ppn=12,pmem=1700M
qsub -I -A wst-164-ab -d `pwd` -l walltime=01:00:00 -q metaq -l nodes=1:ppn=4,pmem=2700m

qsub -l walltime=01:00:00 -l nodes=1:ppn=4,pmem=2700m -N ${fn}_bamstats -q sw -A wst-164-ab

#Access Sequel Pacbio runs via sftp://abacus/lb/robot/pacbioSequel/pacbioRuns
#To log onto abacus type: ssh abacus
#to see modules use: module use ${MUGQIC_INSTALL_HOME}/modulefiles AND THEN module avail AND THEN module load moduleName
http://wiki.genome.mcgill.ca:18080/pages/viewpage.action?pageId=65568820
#If u want to share files with someone u can ask IT to create login credentials for u and they will give u username (eg agiakountis) and password (e.g.  jp98p4wd9).
#Then copy the files to abacus '/lb/scratch/abayega/antonis', login into abacus "ssh abacus" and type
lftp -e "mirror -R  /lb/scratch/abayega/antonis" -u agiakountis sftp://sftp.gqinnovationcentre.com
#Then ask the collaborator to log into gq using the credentials

qsub -A wst-164-ab canu_1-15kb_guillimin.sh
qsub -A wst-164-ab canu_15-30kb_guillimin.sh
qsub -A wst-164-ab canu_40kb-above_guillimin.sh

guillimin and mammoth username = banthony
guillimin and mammoth password = Sip8$@sX
CCRI : vsy-012-01
To sign in: https://portail.calculquebec.ca/accounts/login/
For training http://www.eventbrite.ca/o/calcul-quebec-8295332683
Log on Mammoth:  
ssh banthony@ioannisr-mp2.ccs.usherbrooke.ca

qsub -I -A vsy-012-01 -d `pwd` -l walltime=16:00:00 -q qwork -l nodes=1:ppn=1 #Mammoth mpn=32Gb (24 cores)
qsub -I -A vsy-012-01 -d `pwd` -l walltime=08:00:00 -q qfat256 -l nodes=1:ppn=1 #Mammoth mpn=256Gb (48 cores)

Queue 	Min No of nodes 	Av nodes mpn 	Cores per node 	Maximum run time 	InfiniBand setting
qwork 		1 (24 cores) 	1588 	32 GB 	24 				120 h 					7:2
qfbb 		12 (288 cores) 	216 	32 GB 	24 				120 h 					1:1
qfat256 	1 (48 cores) 	20 		256 GB 	48 				120 h 					1:1
qfat512 	1 (48 cores) 	2 		512 GB 	48 				48 h 					1:1 

#backup_the_mircea_file_again
rsync -av bicyclez@guillimin.hpc.mcgill.ca:/gs/project/wst-164-ab/anthony/Nanopore_data/olive_fly/ /media/abayega/Seagate_U_5TB/analysis/Nanopore_data/olive_fly/

#synchronise!!!! IT WILL DELETE WHATEVER IS ON TARGET FILE BUT NOT ON GUILLIMIN!!!!!!!!!!!
rsync -av --del bicyclez@guillimin.hpc.mcgill.ca:/gs/project/wst-164-ab/spyro /MIRCEA_files_analyse_again/ /lb/project/ioannisr/spyrfs/from_guillimin_march_2017_all_3/MIRCEA_files_analyse_again/

#Error mounting /dev/sdc2 at /media/abayega/New Volume: Command-line `mount -t "exfat" -o "uhelper=udisks2,nodev,nosuid,uid=1000,gid=1000,iocharset=utf8,namecase=0,errors=remount-ro,umask=0077" "/dev/sdc2" "/media/abayega/New Volume"' exited with non-zero exit status 32: mount: unknown filesystem type 'exfat


[bicyclez@sw-2r14-n42 wst-164-aa]$ du -hs anthony/*
3.6T	anthony/Nanopore_data
416K	anthony/Ubuntu progs
12G	anthony/aedes_aegypti
150G	anthony/fruit_fly_illumina_data
32K	anthony/guillimin_submit_clean.sh
3.9G	anthony/hg19
80G	anthony/illumina_reads
128K	anthony/modulefiles
39G	anthony/pacBio_data
6.8G	anthony/software
372M	anthony/spyros_alignment

On guillimin banthony@guillimin.hpc.mcgill.ca is owner 3054995 and group number 13278
On guillimin bicyclez@guillimin.hpc.mcgill.ca is owner 3041357 and group number 13278

uname -a #To determine the type of machine/computer you have in order to choose the right binaries

#On guillimin my $HOME directory is /home/bicyclez
pip3 install --user <package> #To install python packages in my directory without Dunarell you just pip3 install --user <package> and this will install in $HOME/.local/lib/python3.5/site-packages/ or $HOME/local/bin
python2 setup.py install --user
~/.local/bin #python site packages

2to3 #THis is a python tool that converts python2 scripts to python3 scripts

https://broadinstitute.github.io/picard/explain-flags.html #To explain sam flags

perldoc -l XML::Simple #Check is perl module is instaled

pip3 freeze | grep h5py #Find version of installed python module pip freeze lists all installed packages
python -c "import h5py; print h5py.__version__"

#R was misbehaving saying "no such file or directory" and I executed the following to make it work
sudo chown -R $(whoami) ~/.rstudio*

df -h --total #to find space on machine

sed 'NUMq;d' file eg sed '10q;d' 1.txt #to print line 10

tar -cvzf from_my_disk_4th_october_backup.tar.gz ../from_my_disk_4th_october_backup

sed --in-place '20,37d; 45d' file.txt # remove lines 20-37, then line 45 in place

install.packages("readr", lib="/home/banthony/Rpackages") # To install R site packages
library(readr, lib="/home/banthony/Rpackages") # To load the package

install.packages("http://hartleys.github.io/QoRTs/QoRTs_LATEST.tar.gz", repos = NULL, type="source", lib="/home/banthony/Rpackages")

Graham login
ssh banthony@graham.computecanada.ca #Sip8$@sX
salloc --time=24:0:0 --ntasks=2 --account=def-ioannisr

sftp://banthony@graham.computecanada.ca/home/banthony/project/banthony #To access the files
diskusage_report

chown -R banthony:def-ioannisr /home/banthony/projects/def-ioannisr/banthony/analysis/nanopore/olivefly/C004_08A_4_010517/reads/C004_08A_4_010517.tar
chmod a+rwx filex.txt

precompute.sh 1 > ./precompute.000001.out 2>&1


#get all users and their UID
awk -F: '/\/home/ {printf "%s:%s\n",$1,$3}' /etc/passwd

##Docker images
#To see all images
docker images
#To enter an image for example prapi version 1. Use REPOSITORY:TAG
docker run -it prapi:v1 sh
#Once in u can work like in any other terminal
#Mount a folder and operate on it via image
sudo docker run -it --rm -v ${path}:/data prapi:v1


#change file ownership
sudo chown abayega /home/abayega/.Rhistory
sudo chgrp abayega /home/abayega/.Rhistory

#Sharing files with another user on Beluga
https://docs.computecanada.ca/wiki/Sharing_data
In summary, could you please run:
setfacl -m u:"user name":rx "directory to grant access to"
for example:
setfacl -m u:fh4132:rx /scratch/banthony
setfacl -d -m u:fh4132:rx /scratch/banthony/test_d
setfacl -R -m u:fh4132:rX /scratch/banthony/test_d
"-m" means modify permissions, "-d" means to set the default permissions, and "-R" means recursively modify permissions.


