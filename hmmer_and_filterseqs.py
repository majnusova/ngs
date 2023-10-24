## A script that runs HMMER against selected databases and extracts the sequences from the hmmout outfiles (all seqs above the inclusion treshold). Then, it filters the seqs based on their length, X content, and also removes duplicate sequences.

from glob import glob
from Bio import SearchIO
from Bio import SeqIO
import subprocess

# running hhmmsearch against multiple databases at once and saving the resulting files using the .hmmout suffix
def run_hmmsearch(database, hmm_model, output_file):
    command = ["hmmsearch", "-o", output_file, hmm_model, database]
    subprocess.run(command, check=True)

hmm_model = "/home/majnusova/all/projects/bilabrum/data/hmm/MCP_virophages.hmm"
output_directory = "/home/majnusova/all/projects/bilabrum/data/cafeteria/"
database_files = glob("/home/majnusova/all/projects/bilabrum/data/cafeteria/*ORFs_450.fasta") # list of databases 

# Iterate over the list of database files and run hmmsearch 
for database_file in database_files:
    output_file = f"{output_directory}{database_file.split('/')[-1].replace('.fasta', '_MCP.hmmout')}" 
    run_hmmsearch(database_file, hmm_model, output_file) #needs to be nested inside the loop!

# saving IDs of sequences above the inclustion treshold into ids_list
ids_list = []
for file in glob("/home/majnusova/all/projects/bilabrum/data/cafeteria/*MCP.hmmout"):
    hmmer_file = SearchIO.read(file, "hmmer3-text")
    for record in hmmer_file:
        if record.is_included:
            ids_list.append(record.id)
len(ids_list)

# parsing files with ORFs in fasta format
parsed_orfs = {}
for file in glob("/home/majnusova/all/projects/bilabrum/data/cafeteria/*450.fasta"):
    parsed = SeqIO.parse(file, "fasta")
    for orf in parsed:
        parsed_orfs[orf.id] = orf.seq

len(parsed_orfs)

# extrating sequences of our interest from the parsed ORF files 
hmmout_mcp = {}
for orfid, orfseq in parsed_orfs.items():
    if orfid in ids_list:
        hmmout_mcp[orfid] = orfseq
hmmout_mcp

"""
#testing
temp = []
for i in ids_list:
    if i not in hmmout_mcp:
        temp.append(i)
len(temp)
"""

# removing redundancy of the extracted sequences based on the sequence (not ID)
def nonredundant(hmmout_mcp):
    deduplicated_seqs = {}
    temp = []
    for seqid, seq in hmmout_mcp.items():
        if str(seq) not in temp:
            temp.append(str(seq))
            deduplicated_seqs[seqid] = seq
    return deduplicated_seqs
deduplicated = nonredundant(hmmout_mcp)
len(deduplicated)

# filtering based on legth and X proportion 
def filter(deduplicated):
    filtered = {}
    for seqid, seq in deduplicated.items():
        seqlen = len(seq)
        min_len = seqlen
        maxlen = seqlen
        x_count = seq.upper().count("X")
        x_proportion = x_count/seqlen
        if seqlen >= 120 and seqlen <= 850 and x_proportion <= 0.05: # needs to be adjusted for different type of sequences
            filtered[seqid] = str(seq)
    return filtered
finalseqs = filter(deduplicated)
len(finalseqs)

# saving filtered sequences into a file for downstream analyses
with open("/home/majnusova/all/projects/bilabrum/results/mcp_filtered.fas", "w") as outfile:
    for seqid, sequence in finalseqs.items():
        outfile.write(f">{seqid}\n{sequence}\n")

