#!/bin/sh
#PBS -N index_genome
#PBS -q batch 
#PBS -l nodes=1:ppn=80 
#PBS -l walltime=999:99:99 

source activate snakemake

hisat2-build /home/users/mvladka/projects/jotnarlogs/data/aurantiochytrium/Aurli1_AssemblyScaffolds.fasta /home/users/mvladka/projects/jotnarlogs/intermediate_files/aurli_index &> /home/users/mvladka/projects/jotnarlogs/scripts/aurli_indexing.log
