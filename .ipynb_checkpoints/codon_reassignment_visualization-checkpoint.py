from Bio.Blast import NCBIXML
from Bio import SeqIO

translation_table = {"TTT": "F", "TTC": "F", "TTA": "L" , "TTG": "L",
                "TCT": "S", "TCC": "S" , "TCA": "S" , "TCG": "S",
                "TAT": "Y" , "TAC": "Y" , "TAA": "*" , "TAG": "q",
                "TGT": "C" , "TGC": "C" , "TGA": "*" , "TGG": "W",
                "CTT": "L" , "CTC": "L" , "CTA": "L" , "CTG": "L",
                "CCT": "P" , "CCC": "P" , "CCA": "P" , "CCG": "P",
                "CAT": "H" , "CAC": "H" , "CAA": "Q" , "CAG": "Q",
                "CGT": "R" , "CGC": "R" , "CGA": "R" , "CGG": "R",
                "ATT": "I" , "ATC": "I" , "ATA": "I" , "ATG": "M",
                "ACT": "T" , "ACC": "T" , "ACA": "T" , "ACG": "T",
                "AAT": "N" , "AAC": "N" , "AAA": "K" , "AAG": "K",
                "AGT": "S" , "AGC": "S" , "AGA": "R" , "AGG": "R",
                "GTT": "V" , "GTC": "V" , "GTA": "V" , "GTG": "V",
                "GCT": "A" , "GCC": "A" , "GCA": "A" , "GCG": "A",
                "GAT": "D" , "GAC": "D" , "GAA": "E" , "GAG": "E",
                "GGT": "G" , "GGC": "G" , "GGA": "G" , "GGG": "G"}

blastout = open('/home/vmaj/programs/blast/blastout_atp.txt') #tblastn, proteinove query = ATPazy z bilabra, database = concatenated asseblies -> parsujeme tento blast output
blast_records = NCBIXML.parse(blastout)
x = next(blast_records)

x.query

nucl_dict = {} #ids a nukl sekvence vsech nodu ze vsech assemblies
for record in SeqIO.parse("/home/vmaj/programs/blast/combined_assembly.fas", "fasta"):
    nucl_dict[record.id] = record.seq
nucl_dict

def rev_com(seq):
    seq = seq.upper()
    seq = seq.replace('A', 't')
    seq = seq.replace('T', 'a')
    seq = seq.replace('G', 'c')
    seq = seq.replace('C', 'g')
    result = seq.upper()
    return result[::-1]

id_nucl = {} #ids + nukl sekvence proteinu, ktere nas zajimaji
blastout = open('/home/vmaj/programs/blast/blastout_atp.txt')
blast_records = NCBIXML.parse(blastout)
for record in blast_records:
    query = record.query
    best_blast_res = record.alignments[0].hsps[0]
    hit_id = record.alignments[0].hit_id
    start = best_blast_res.sbjct_start
    end = best_blast_res.sbjct_end
    frame = best_blast_res.frame[1]
    if frame > 0:
        id_nucl[query] = nucl_dict[hit_id][start-1:end+3]
    else:
        rev_start = start-4
        if rev_start < 0:
            rev_start = start - 1
        seq = nucl_dict[hit_id][rev_start:end]
        rev_comp_seq = rev_com(seq)
        id_nucl[query] = rev_comp_seq
id_nucl

with open('/home/vmaj/programs/hmmer3/cds.txt', 'w') as f:
    for node, cds in id_nucl.items():
        f.write(f'>{node}\n{cds}\n')
        
translated_seqs = {}

for id_node, nucl_seq in id_nucl.items():#id_nucl osahuje ID:DNA seq
    protein = ""
    for nt in range(0, len(nucl_seq)-2, 3):
        codon = nucl_seq[nt : nt+3]
        protein += translation_table[codon]
        translated_seqs[id_node] = protein
        #return id_node
translated_seqs


with open('/home/vmaj/programs/hmmer3/translated_seqs_q.txt', "w") as outfile:
    for id_node, prot_seq in translated_seqs.items():
        outfile.write(f">{id_node}\n{prot_seq}\n")

