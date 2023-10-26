#!/bin/sh
#PBS -N snakejob
#PBS -q batch
#PBS -l nodes=1:ppn=80
#PBS -l walltime=999:00:00

cd /home/users/mvladka/projects/jotnarlogs/cilia_snake

source activate snakemake

snakemake --unlock -s workflows/Snakefile
snakemake --cores 80 -s workflows/Snakefile
