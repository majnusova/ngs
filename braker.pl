# How to run braker3 on genome.osu.cz:

# if not done yet, you need to copy augustus config file to a writable location and export its path:
# cp -r /opt/miniforge3/envs/braker3/config/ /home/mvladka/annotations
# nano ~/.bashrc and add this variable: export AUGUSTUS_CONFIG_PATH=/home/mvladka/annotations
# source ~/.bashrc
# test if the path was exported correctly: echo $AUGUSTUS_CONFIG_PATH (it should print this path: /home/mvladka/annotations)

# -------------------------
screen -S braker
source ~/.bashrc # it probably has to be done in everytime you start a new screen session
conda activate braker3 # activate the env within the screen
braker.pl --genome=/home/mvladka/annotations/asm_k77_raw.gc_lessN.scafSeq.fasta --species=pseudellepsioidon --rnaseq_sets_ids=pseuedell-edaph_trim --rnaseq_sets_dirs=/home/mvladka/annotations --prot_seq=/home/mvladka/annotations/Stramenopiles.fa --threads=8 &> anot.log 

# --rnaseq_sets_ids=pseuedell-edaph_trim needs to be set like this if you have paired-end reads, for unapaired reads: --rnaseq_sets_ids=id1,id2

#threads setting:
"""
the total number of threads used by your jobs should never exceed 16. As for 10 threads - that is usually quite fine.
You can use 'top'/'htop' to see overall system load."""
# if you have a highly fragmented genome, only one thread should be used (linear mode)
