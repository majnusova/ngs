from Bio import SeqIO
import subprocess
import os
import glob
#subprocess.run('source activate hhsuite', shell=True)

# parsing mcp.fas which contains sequences that serve as inputs for hhblits
def parse_seqs():
    with open("/home/users/mvladka/hhpred_databases/mcp.fas", "r") as infile:
        parsed_seqs = {}
        for seq in SeqIO.parse(infile, "fasta"):
            parsed_seqs[seq.id] = seq.seq #str()
    return parsed_seqs
input_seqs = parse_seqs()

# saving each sequence into a separate temporary file (hhblits accepts only files)
outdir = "/home/users/mvladka/hhpred_databases/inputs_results"
for seqid, seq in input_seqs.items():
    filename = os.path.join(outdir, f"{seqid}.fa")
    with open(filename, "w") as outfile:
        outfile.write(f">{seqid}\n{seq}\n")

# defining the hhblits command (these 4 databases are useful for divergent viral seqs)
def hhblits(input_file, output_file):
    command = [
        "hhblits",
        "-i", input_file,
        "-d", "COG_KOG",
        "-d", "pdb70",
        "-d", "NCBI_CD",
        "-d", "pfam",
        "-o", output_file,
        "-cpu", "8",
        "-n", "3"
    ]
    subprocess.run(command, check=True)

# Iterate over the input seq files and run hhblits with all the sequences one by one
for filename in glob.glob("/home/users/mvladka/hhpred_databases/inputs_results/*.fa"):
    #output_directory = "/home/users/mvladka/hhpred_databases/results_hhblits"
    input_file = filename
    output_file = f"{filename}.hhr"
    hhblits(input_file, output_file)

# removing the temporary .fa files we created
for filename in glob.glob("/home/users/mvladka/hhpred_databases/inputs_results/*.fa"):
    os.remove(filename)
