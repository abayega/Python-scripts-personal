#!/bin/bash
#SBATCH --account=def-ioannisr
#SBATCH --cpus-per-task=24      # number of MPI processes
#SBATCH --time=0-12:00          # time (DD-HH:MM)
#SBATCH --job-name=tar
#SBATCH --mem=40G

cd /home/banthony/scratch/reads/pacbio/olivefly/dnaseq


#awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' olivefly_S1_L002_R1_001.trim16.filtpaired.fastq > olivefly_S1_L002_R1_001.trim16.filtpaired.fasta && awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' olivefly_S1_L002_R2_001.filtpaired.fastq > olivefly_S1_L002_R2_001.filtpaired.fasta

tar -zcvf boleae_pacbio_male.merged.filtered_subreads.tar.gz boleae_pacbio_male.merged.filtered_subreads
