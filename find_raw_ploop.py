from glob import glob
import re
import subprocess
from Bio import SearchIO
from Bio import SeqIO
from Bio.Blast import NCBIXML


def run_hmmsearch(query_file, hmm_database, output_file):
    # Build the command
    command = ["hmmsearch", "-o", output_file, hmm_database, query_file]
    # Execute the command
    subprocess.run(command, check=True)


# ### Creating a dictionary: keys are the names of the organisms to check for RAW and values are paths to the Eukprot files of these orgs.

eukprot_raw = glob("/home/vmaj/Eukprot_sequences/*.fasta")
eukprot_dict = dict() # {'Muranothrix_gubernata': '/home/vmaj/Eukprot_sequences/EP00824_Muranothrix_gubernata.fasta',...}

for file_name in eukprot_raw:
    name = file_name.split("_")[2:]
    name = "_".join(name)[:-6]
    eukprot_dict[name] = file_name

len(eukprot_dict)


# ### Running hmmsearch: raw_hmm against selected eukprot files and saving the -hmmout results.

targets = open('/home/vmaj/PROJECTS/Jotnarlogs/distribution_mapping/orgs_to_map_names.txt')
hmm_database = "/home/vmaj/PROJECTS/Jotnarlogs/distribution_mapping/raw_hmm.txt"

for target in targets:
    target = target.strip()
    if target in eukprot_dict:
        out_folder = "/home/vmaj/PROJECTS/Jotnarlogs/distribution_mapping/hmm_results/"
        out_file = f"{out_folder}{target}.hmmout"
        query_file = eukprot_dict[target]
        run_hmmsearch(query_file, hmm_database, out_file) #running hmmsearch


# ### Parsing hmmer outputs and filtering only those results that have e-value above e-10. hsp midline??

results_raw = glob('/home/vmaj/PROJECTS/Jotnarlogs/distribution_mapping/hmm_results/*.hmmout')
best_hits_aln = {}
for res_file in results_raw:
    res = SearchIO.parse(res_file, 'hmmer3-text')
    for hit in res:
        for hsp in hit.hsps:
            evalue_exponent = float(str(hsp.evalue).split("e")[-1]) #make string from evalue, split it on "e" to take the exponent and turn it back to float
            if evalue_exponent <= -10:
                query_aligned = hsp.aln[0].seq  #Access query aligned sequence
                hit_aligned = hsp.aln[1].seq  #Access hit aligned sequence
                best_hits_aln[hsp.hit.id] = {'query_aligned': query_aligned, 'hit_aligned': hit_aligned}

best_hits_aln # {'EP00023_Dictyostelium_discoideum_P009877': {'query_aligned': Seq('asaaqrvvkvavlGdarsGkaslvrrlitaefdeqyletlgievselaaleads...yvl'),'hit_aligned': Seq('KFDNTKEIKLVLIGDGGVGKSTYINRLLTGEFETQYVATFGCSVHKFNFKTTIG...KLI')},

with open('/home/vmaj/PROJECTS/Jotnarlogs/distribution_mapping/wtf.txt', 'w') as outfl:
    for i in best_hits_aln:
        outfl.write(f'{i}\n')
# ### Searching for p-loop motif GxxxxG in the aligned region of the HIT and QUERY - maybe possible to use hsp_midline? zkusit udelat zpetny blast rovnou bez techto podminek pro kontrolu

"""checks if the match object returned by re.search() is not None. In Python, non-empty objects are considered True, and None is considered False. 
So if match is equivalent to if match is not None.
The if match in syntax is not valid in this context because match is not a collection or iterable object that can be used with the in operator. 
The re.search() function returns a single match object or None, and we want to check if a match is found (match is not None) rather than checking if match exists in a specific collection."""
putative_raz = {}
p_loop = r'G.{4}GK'
for hit_id, aligned_seqs in best_hits_aln.items():
    query_aligned = str(aligned_seqs['query_aligned']).upper()
    hit_aligned = str(aligned_seqs['hit_aligned']).upper()

    hit_match = re.search(p_loop, hit_aligned) #searches for p_loop in the aligned region of the hit
    query_match = re.search(p_loop, query_aligned)

    if hit_match and query_match: # if hit_macth and query_match not None
        ploop_hit_start = hit_match.start() # searches for the start of the ploop in the hit seq
        ploop_hit_end = hit_match.end()
        ploop_query_start = query_match.start()
        ploop_query_end = query_match.end()

        # coordinates of the p-loop motif
        p_loop_hit_aligned = hit_aligned[ploop_hit_start:ploop_hit_end]
        p_loop_query_aligned = query_aligned[ploop_query_start:ploop_query_end]

        """nezastavi se to u prvniho C?"""
        # Check if any "C" occurs after the p_loop motif in both hit_aligned and query_aligned. If so, an index is added to cysteines and the indices are appended to the lists
        hit_positions_c = []
        query_positions_c = []
        for index, aa in enumerate(hit_aligned[ploop_hit_end:]):
            if aa == "C":
                hit_positions_c.append(index + ploop_hit_end) #getting the exact position of the C
        
        for index, aa in enumerate(query_aligned[ploop_query_end:]):
            if aa == "C":
                query_positions_c.append(index + ploop_query_end)

        # Check if any of the "C" occurrences are at the same position in hit_aligned and query_aligned
        found_aligned_c = False # to initialize the variable found_common_c, boolean variable False
        for pos in hit_positions_c:
            if pos in query_positions_c:
                found_aligned_c = True #variable is se to True if found
                # break # break statement is executed if shared position found, it terminates the loop, skipping the remaining iterations - > if found_alignd_c: is executed

        if found_aligned_c:
            putative_raz[hit_id] = {'query_aligned': query_aligned, 'hit_aligned': hit_aligned}

putative_raz



# ### Extracting IDs of putative RAZ sequences and getting the names of the respective eukprot files.

#IDs of protein seqs
eukprot_ids = []
for euk_id in putative_raz.keys():
    eukprot_ids.append(euk_id)
eukprot_ids

#corresponding file names - toto uz je udelane v eukprot_dict
eukprot_files = set()
for euk_id in putative_raz.keys():
    general_name = euk_id.split("_")[1:]
    general_name = general_name [:-1]
    general_name = "_".join(general_name)
    eukprot_files.add(general_name)
eukprot_files


# ### Parsing Eukprot files of our interest. prepsat, uz je v eukprot:dict


paths = []
for euk_file in eukprot_files:
    euk_file = euk_file.strip() 
    if euk_file in eukprot_dict:
        paths.append(eukprot_dict[euk_file])
paths

parsed_eukfiles = {}
for path in paths:
    for record in SeqIO.parse(path, "fasta"):
        descr = record.id
        seq = record.seq
        parsed_eukfiles[descr] = seq
len(parsed_eukfiles)


# ### Extracting putative RAZ sequences from eukprot to be able to use them for control blast.


raz_candidates = {}
for prot_id in eukprot_ids:
    for descr, seq in parsed_eukfiles.items():
        if prot_id in descr:
            raz_candidates[descr] = seq
raz_candidates    


# ### Saving putative RAZ sequences (Eukprot ID + seq) into a file that would serve as a query file for blastp.

with open("/home/vmaj/PROJECTS/Jotnarlogs/distribution_mapping/control_blast_results/raz_candidates.fas", "w") as outf:
    for i, j in raz_candidates.items():
        outf.write(f">{i}\n{j}\n")


# ### Running blastp for humans: candidate RAZ seqs x rabs database.

#do klasickeho text file ukladam vysledky zpetneho blast proti rabs database
database = "/home/vmaj/programs/blast/rabs_final.txt"
query_file = "/home/vmaj/PROJECTS/Jotnarlogs/distribution_mapping/control_blast_results/raz_candidates.fas"
output_dir = "/home/vmaj/PROJECTS/Jotnarlogs/distribution_mapping/control_blast_results"

cmd = ["blastp", "-query", query_file, "-db", database, "-out", f"{output_dir}/control_blast_results.txt"]
subprocess.run(cmd, check=True, text=True)


# ### Running blastp for computer: candiate RAZ seqs x rabs database.

database = "/home/vmaj/programs/blast/rabs_final.txt"
query_file = "/home/vmaj/PROJECTS/Jotnarlogs/distribution_mapping/control_blast_results/raz_candidates.fas"
output_dir = "/home/vmaj/PROJECTS/Jotnarlogs/distribution_mapping/control_blast_results"

output_file = f"{output_dir}/control_blast_results.xml"
cmd = ["blastp", "-query", query_file, "-db", database, "-outfmt", "5", "-out", output_file] #xml file
subprocess.run(cmd, check=True, text=True)


# ### Saving IDs of RAZ sequences into a list.

blast_result_file = "/home/vmaj/PROJECTS/Jotnarlogs/distribution_mapping/control_blast_results/control_blast_results.xml"

raw_sequences = set()

blast_records = NCBIXML.parse(open(blast_result_file))
for blast_record in blast_records:
    query_sequence = blast_record.query

    alignment_count = 0  #counter for the number of alignments processed (nekdy je raw az uplne dole a jen jednou - nebude to ono)

    for alignment in blast_record.alignments:
        alignment_count += 1
        if alignment_count > 20:  # Check if more than 20 alignments have been processed
            break

        for hsp in alignment.hsps:
            hit_def = alignment.title
            #if re.findall(r'raw', hit_def, re.IGNORECASE):
            if "raw" in hit_def.lower():
                raw_sequences.add(query_sequence)
                break  # Skip remaining alignments for the current query sequence

raw_sequences


#uolozit si raw_sequences i se seqs do file

{'EP00002_Diphylleia_rotans_P023087', ANO
 'EP00002_Diphylleia_rotans_P023088', ANO
 'EP00033_Pygsuia_biforma_P002446', ANO
 'EP00113_Trichoplax_adhaerens_P006561', NE
 'EP00130_Batrachochytrium_dendrobatidis_P006821', ANO
 'EP00131_Spizellomyces_punctatus_P005997', ANO
 'EP00157_Rozella_allomycis_P002373', ANO
 'EP00473_Plasmodiophora_brassicae_P006349', ANO
 'EP00705_Trichomonas_vaginalis_P030997', SPIS NE
 'EP00735_Rhodelphis_limneticus_P004029', ANO
 'EP00736_Rhodelphis_marinus_P001758', ANO
 'EP00741_Cyanophora_paradoxa_P010914', ANO
 'EP01135_Chromosphaera_perkinsii_P007894', ANO
 'EP01137_Pigoraptor_chileana_P018232', ANO
 'EP01137_Pigoraptor_chileana_P018233'} ANO
# ### Extract RAZ protein sequences from Eukprot files:
raw_sequences = ['EP00002_Diphylleia_rotans_P023087',
 'EP00002_Diphylleia_rotans_P023088',
 'EP00033_Pygsuia_biforma_P002446',
 'EP00113_Trichoplax_adhaerens_P006561',
 'EP00130_Batrachochytrium_dendrobatidis_P006821',
 'EP00131_Spizellomyces_punctatus_P005997',
 'EP00157_Rozella_allomycis_P002373',
 'EP00473_Plasmodiophora_brassicae_P006349',
 'EP00705_Trichomonas_vaginalis_P030997',
 'EP00735_Rhodelphis_limneticus_P004029',
 'EP00736_Rhodelphis_marinus_P001758',
 'EP00741_Cyanophora_paradoxa_P010914',
 'EP01135_Chromosphaera_perkinsii_P007894',
 'EP01137_Pigoraptor_chileana_P018232',
 'EP01137_Pigoraptor_chileana_P018233']
raw_seqs_extracted = {}
for prot in raw_sequences:
    for descr, seq in parsed_eukfiles.items():
        if prot in descr:
            raw_seqs_extracted[descr] = seq
    
   
