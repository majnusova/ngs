#!/bin/sh
#PBS -N hhblits_subprocess
#PBS -q batch
#PBS -l nodes=1:ppn=80
#PBS -l walltime=999:00:00

source activate hhsuite
python3 /home/users/mvladka/hhpred_databases/runhhblits.py &> hhblits.log
