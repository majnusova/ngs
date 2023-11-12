#!/bin/sh
#PBS -N virsorter2
#PBS -q batch 
#PBS -l nodes=1:ppn=80 
#PBS -l walltime=999:99:99 
source activate /home/users/mvladka/miniconda3/envs/vs2
			
cd /home/users/mvladka/projects/virsorter

virsorter run -w bilabrum_virsorter2.out -i bilabrum_metaspades_miseq.fasta --include-groups "dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae" --min-length 450 -j 4 &> vs2_bil.log
