from glob import glob
from Bio import SearchIO
from Bio import SeqIO
import subprocess

# sequences from the dataset by bellas and sommaruga contained spaces
def fix_fasta_headers(input_file, output_file):
    with open(output_file, "w") as outfile:
        for record in SeqIO.parse(input_file, "fasta"):
            new_id = record.description.replace(" ", "_")
            outfile.write(f">{new_id}\n{record.seq}\n")

input_file = "/home/majnusova/all/projects/bilabrum/results/clans/mcp_final_seqs_for_cyto.fasta"
output_file =  "/home/majnusova/all/projects/bilabrum/results/clans/mcp_final_seqs_for_cyto_fixed.fasta"

fix_fasta_headers(input_file, output_file)


file = open("/home/majnusova/all/projects/bilabrum/results/clans/mcp_final_seqs_for_cyto_fixed.fasta")
parsed_file = SeqIO.parse(file, "fasta")
def x_proportion(parsed_file):
    x_proportions = {}
    for protein in parsed_file:
        xcount = protein.seq.upper().count("X")
        seq_length = len(protein.seq)
        proportion = (xcount/seq_length)
        if proportion <= 0.05:
            x_proportions[protein.id] = protein.seq
    return x_proportions
x_proportion_in_seq = x_proportion(parsed_file)
print(len(x_proportion_in_seq))


def remove_partial_seqs(x_proportions):
    dict_150aa = {} # id : seqs 
    for prot_id, prot_sequence in x_proportions.items():
        min_length = len(prot_sequence)
        max_length = len(prot_sequence)
        if min_length >= 130 and max_length <= 1000: # nejdelsi pATPasa chteneho typu ma 344aa; nejdelsi MCP 614 aa
            dict_150aa[prot_id] = prot_sequence
    return dict_150aa
res = remove_partial_seqs(x_proportion_in_seq)
print(len(res))


def remove_redundant_seqs(res):
    temp = []
    deduplicated_seqs = {}
    for prot_name, prot_seq in res.items():
        if prot_seq not in temp:
            temp.append(prot_seq)
            deduplicated_seqs[prot_name] = prot_seq
    return deduplicated_seqs
res2 = remove_redundant_seqs(res)
print(len(res2))


# saving deduplicated and filtered sequences into a file
def save_results(res2): 
    with open("/home/majnusova/all/projects/bilabrum/results/clans/deduplicated_MCP_for_cytoscape_final.txt", "w") as outfile:
        for prot_name, prot_seq in res2.items():
            outfile.write(f'>{prot_name}\n{prot_seq}\n')
res3 = save_results(res2)        

