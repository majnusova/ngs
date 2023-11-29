#!/bin/sh
#PBS -N trim
#PBS -q batch 
#PBS -l nodes=1:ppn=80 
#PBS -l walltime=999:99:99 


source activate /home/users/mvladka/miniconda3/envs/assembly-env
			
cd /home/users/mvladka/projects/pseudellipsoidon

trimmomatic PE -threads 20 pseuedell-edaph_1.fastq pseuedell-edaph_2.fastq -baseout pseuedell-edaph ILLUMINACLIP:/home/users/mvladka/miniconda3/envs/assembly-env/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 &> pseud_trim.log
