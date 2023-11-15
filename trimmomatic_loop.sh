#!/bin/sh
#PBS -N trim
#PBS -q batch 
#PBS -l nodes=1:ppn=80 
#PBS -l walltime=999:99:99 
#PBS -j oe

source activate /home/users/mvladka/miniconda3/envs/assembly-env
			
cd /home/users/mvladka/projects/assemblies/bilabrum


# Trimmomatic parameters
ADAPTERS="/home/users/mvladka/miniconda3/envs/assembly-env/share/trimmomatic/adapters/TruSeq3-PE-2.fa"
LEADING="3"
TRAILING="3"
SLIDINGWINDOW="4:15"
MINLEN="75"

# Loop through each pair of files
for file in *_1.fastq.gz
do
    # Names for input files
    file1=${file}
    file2=${file/_1.fastq.gz/_2.fastq.gz}

    # Base name for output
    base=$(basename ${file} _1.fastq.gz)

    # Run Trimmomatic
    trimmomatic PE -threads 30 \
        ${file1} ${file2} \
        ${base}_1_paired.fastq.gz ${base}_1_unpaired.fastq.gz \
        ${base}_2_paired.fastq.gz ${base}_2_unpaired.fastq.gz \
        ILLUMINACLIP:${ADAPTERS}:2:20:10 \
        LEADING:${LEADING} \
        TRAILING:${TRAILING} \
        SLIDINGWINDOW:${SLIDINGWINDOW} \
        MINLEN:${MINLEN} &> ${base}_trim_log
done

