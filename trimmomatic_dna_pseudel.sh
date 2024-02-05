#!/bin/sh
#PBS -N trim
#PBS -q batch
#PBS -l nodes=1:ppn=80
#PBS -l walltime=999:99:99


source activate /home/users/mvladka/miniconda3/envs/assembly-env

cd /home/users/mvladka/projects

trimmomatic PE -threads 20 lane114s005734_1_sequence.fastq lane114s005734_2_sequence.fastq -baseout lane114s005734 ILLUMINACLIP:/home/users/mvladka/miniconda3/envs/assembly-env/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:20:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 &> pseud_dna_trim.log

"""
#TruSeq3-PE-2.fa is one of the adapter sequence files provided with Trimmomatic
# 2:20:10 =  2: Seed Mismatches. This number allows up to 2 mismatches between the read sequence and the adapter sequence; allowing mismatches helps to identify and remove adapters even when there are small discrepancies due to sequencing errors or variations
# 20: Palindrome Clip Threshold. This threshold is used for detecting adapter sequences that are inserted as palindromes, a situation that can occur when very short fragments are sequenced. The threshold value of 20 specifies the minimum alignment score needed to clip the adapter sequence. The alignment score is calculated based on matches between the read and the adapter sequence, with mismatches penalized. A higher score means a higher confidence in the adapter sequence's presence.
# 10: Simple Clip Threshold. This value is used for straightforward (non-palindromic) adapter clipping. If the alignment score between the adapter sequence and any part of the read exceeds 10, the adapter sequence is clipped off. This score is simpler than the palindrome clip threshold, focusing on direct matches between the adapter and read sequences.

#LEADING:3: Removes leading low-quality bases (below quality 3).
#TRAILING:3: Removes trailing low-quality bases (below quality 3).
#SLIDINGWINDOW:4:15: Scans the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15.
#MINLEN:75: Drops reads which are less than 75 bases long after all the above trimming operations.
"""
