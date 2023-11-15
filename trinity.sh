#!/bin/sh
#PBS -N caecitellus
#PBS -q batch 
#PBS -l nodes=1:ppn=80 
#PBS -l walltime=999:99:99 

source activate /home/users/mvladka/miniconda3/envs/assembly-env
			
cd /home/users/mvladka/projects/assemblies/bilabrum

Trinity --trimmomatic --seqType fq --max_memory 120G --left cae1_SRR24392492_1_paired.fastq.gz,cae2_SRR24392493_1_paired.fastq.gz,cae3_SRR24392494_1_paired.fastq.gz --right cae1_SRR24392492_2_paired.fastq.gz,cae2_SRR24392493_2_paired.fastq.gz,cae3_SRR24392494_2_paired.fastq.gz --CPU 75 --full_cleanup --output trinity_out &> trinity_errors.log



