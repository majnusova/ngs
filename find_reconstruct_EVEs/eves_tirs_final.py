#!/usr/bin/env python
# coding: utf-8

# In[62]:


#!/usr/bin/env python
# coding: utf-8

# ### to run this code, you need to have standalone EMBOSS (getorf), BLAST, and HMMER installed



import subprocess
import os
from glob import glob
from Bio import SearchIO
from Bio import SeqIO 
from Bio.Blast import NCBIXML 
from collections import defaultdict


# ## load and adjust all paths
# /home/majnusova/all/projects/plv/paratrimastix/ - vzdy se jen nahradi tato cast cesty v celem skriptu

# Set absolute path to your hmm model and genomes - these locations need to exist already!:
hmm_model = "/home/majnusova/all/projects/plv/paratrimastix/polb_classI.hmm"
input_genomes = glob("/home/majnusova/all/projects/plv/paratrimastix/input_genomes/*.fasta") #genomes need to end in .fasta suffix



# ### Extract ORFs from genomes
def run_getorf(genome, orf_file):
    command = ["getorf", "-sequence", genome, "-outseq", orf_file, "-table", "0", "-minsize", "450", "-find", "1"] 
    subprocess.run(command, check=True)

outdir = "/home/majnusova/all/projects/plv/paratrimastix/getorf_genomes/" 
if not os.path.exists(outdir):
    os.makedirs(outdir)

for genome in input_genomes:
    orf_file = f"{outdir}{genome.split('/')[-1].replace('.fasta', '_orfs.fa')}"
    run_getorf(genome, orf_file)


# In[63]:


# ### hmmsearch with HMM against those extracted ORFs -> hmmout
# running hhmmsearch against multiple databases at once and saving the resulting files using the .hmmout suffix
def run_hmmsearch(orf_file, hmm_model, output_file):
    command = ["hmmsearch", "-o", output_file, hmm_model, orf_file]
    subprocess.run(command, check=True)

output_directory = "/home/majnusova/all/projects/plv/paratrimastix/hmmout/" 
if not os.path.exists(output_directory):
    os.makedirs(output_directory)
orf_files = glob(os.path.join(outdir, "*_orfs.fa"))


# iterate over the list of database files and run hmmsearch 
for orf_file in orf_files:
    output_file = f"{output_directory}{orf_file.split('/')[-1].replace('_orfs.fa', '.hmmout')}" 
    run_hmmsearch(orf_file, hmm_model, output_file) #needs to be nested inside the loop!


# In[64]:


# ### saving IDs of ALL homologs above inclusion threshold found by HMMER and name of the assembly they were found in
# ### genome_scaffold_ids - genome: ids
# saving IDs of sequences above the inclustion treshold into ids_list
ids_list_orfs = []
for file in glob("/home/majnusova/all/projects/plv/paratrimastix/hmmout/*.hmmout"):
    hmmer_file = SearchIO.read(file, "hmmer3-text")
    for record in hmmer_file:
        if record.is_included:
            ids_list_orfs.append(record.id)
ids_list_orfs
# In[ ]:


# In[65]:


# names of viral orfs and the genome they originate from
genome_scaffold_orfs = {}

for file in glob("/home/majnusova/all/projects/plv/paratrimastix/hmmout/*.hmmout"):
    genome_name = os.path.basename(file).split('.hmmout')[0] # getting the real names of the genomes; os.path.basename(file) extracts the base name from the file's full path.
    if genome_name not in genome_scaffold_orfs:
        genome_scaffold_orfs[genome_name] = []

    hmmer_file = SearchIO.read(file, "hmmer3-text")
    
    for record in hmmer_file:
        if record.is_included:
            scaffold_id = record.id  # scaffold ID - removing _orf
            genome_scaffold_orfs[genome_name].append(scaffold_id)
            
total_viral_orfs = sum(len(i)for i in genome_scaffold_orfs.values())
genome_scaffold_orfs 


# In[66]:


# creating a dict of lists (genomes: ids)
genome_scaffold_ids = {}

for file in glob("/home/majnusova/all/projects/plv/paratrimastix/hmmout/*.hmmout"):
    genome_name = os.path.basename(file).split('.hmmout')[0] # getting the real names of the genomes; os.path.basename(file) extracts the base name from the file's full path.
    if genome_name not in genome_scaffold_ids:
        genome_scaffold_ids[genome_name] = []

    hmmer_file = SearchIO.read(file, "hmmer3-text")
    
    for record in hmmer_file:
        if record.is_included:
            scaffold_id = record.id.rsplit("_", 1)[0] # scaffold ID - removing _orf
            genome_scaffold_ids[genome_name].append(scaffold_id)



total_ids = sum(len(id) for id in genome_scaffold_ids.values())


# ### Extracting scaffolds with the viral homologs found by HMMER (ids_list), scaffolds shorter than 18000 nt are discarded
# ids_list contains IDs of ORFs, not scaffolds! this list needs to be modified accordingly
ids_list_scaffolds = []
for id in ids_list_orfs:
    id_modified = id.rsplit('_', 1)[0] # removing the number of orf -> id of the scaffold
    ids_list_scaffolds.append(id_modified)
#print(ids_list_scaffolds)
# #### ids_list_scaffolds = all viral scaffolds above inclusion threshold (may be shorter than 20000)
#print(ids_list_scaffolds)
# ## to be able to extract the scaffolds, blastable databases need to be created first (for genome and orfs)
genome_scaffold_ids


# In[67]:


# In[ ]:


def makeblastdb(file_path, db_type, blastables_outdir):
    db_name = f"{blastables_outdir}/{file_path.split('/')[-1].split('.')[0]}"
    command = ["makeblastdb", "-in", file_path, "-dbtype", db_type, "-out", db_name, "-parse_seqids"]
    subprocess.run(command, check=True)

blastables_outdir = "/home/majnusova/all/projects/plv/paratrimastix/blastables/"
if not os.path.exists(blastables_outdir):
    os.makedirs(blastables_outdir)


for genome_path in input_genomes:
    makeblastdb(genome_path, "nucl", blastables_outdir) # misto file_path volame s genome_path

for orf_file_path in orf_files:
    makeblastdb(orf_file_path, "prot", blastables_outdir)


# In[68]:


# ### modifying a dictionary for blastdbcmd
# #### genome_scaffold_orfs_descr contains the same IDs as genome_scaffold_ids + description


# creating a dict of lists (genomes: ids)
# this dict contains also the description - it allows the usage of the -range option of blastdbcmd command
genome_scaffold_orfs_descr = {}

for file in glob("/home/majnusova/all/projects/plv/paratrimastix/hmmout/*.hmmout"):
    genome_name = os.path.basename(file).split('.hmmout')[0] # getting the real names of the genomes; os.path.basename(file) extracts the base name from the file's full path.
    if genome_name not in genome_scaffold_orfs_descr:
        genome_scaffold_orfs_descr[genome_name] = []

    hmmer_file = SearchIO.read(file, "hmmer3-text")
    
    for record in hmmer_file:
        if record.is_included:
            genome_scaffold_orfs_descr[genome_name].append(record.id + ";" + record.description)



genome_scaffold_orfs_descr


# In[69]:


total_ids = sum(len(id) for id in genome_scaffold_orfs_descr.values())
total_ids


# In[70]:


# ### Extracting scaffolds (30000nt downstream and upstream of the viral genes) containing viral genes detected by HMMER. 
# ### Scaffolds shorter than 20000 nts are discarded.



scaffolds_of_interest = []

def run_blastdbcmd(id, genome_db, outscaffold, range):
    command = ["blastdbcmd", "-entry", str(id), "-db", genome_db, "-out", outscaffold, "-range", range]
    print(f"Running command: {' '.join(command)}")  

    result = subprocess.run(command, text=True, capture_output=True)

    if result.returncode != 0:
        print(f"Warning: blastdbcmd failed for ID {id} in database {genome_db}. Error: {result.stderr}")
        if os.path.exists(outscaffold):
            os.remove(outscaffold)
        return False, None

    if os.path.exists(outscaffold) and os.path.getsize(outscaffold) > 0:
        record = SeqIO.read(outscaffold, "fasta")
        if len(record.seq) < 20000:
            print(f"Notice: Scaffold {id} is shorter than 20000 nucleotides, removing file.")
            os.remove(outscaffold)
            return False, None
        else:
            scaffolds_of_interest.append(id)
            return True, len(record.seq)
    else:
        print(f"Notice: No scaffold file created for ID {id}.")
        return False, None

scaffolds_dir = "/home/majnusova/all/projects/plv/paratrimastix/eve_outscaffolds/"
if not os.path.exists(scaffolds_dir):
    os.makedirs(scaffolds_dir)
    
blastables_outdir = "/home/majnusova/all/projects/plv/paratrimastix/blastables/"
prot_database = "_orfs."

for genome, ids in genome_scaffold_orfs_descr.items(): #ids = NODE_15445_length_2031_cov_23.9084_2;[1007 - 1726] No definition line found
    print(f"Processing genome: {genome}")  
    genome_dir = os.path.join(scaffolds_dir, genome)
    os.makedirs(genome_dir, exist_ok=True)
    for i in ids:
        res_scaf_orf = i.split(";")[0]
        id_full, descr = i.split(';')
        id = id_full.rsplit('_', 1)[0]  
        print(f"Processing ID: {id}")  

        if "REVERSE" in descr:
            start_orf_part = descr.split(" - ")[1]
            start_orf = start_orf_part.split("]")[0]
            end_orf_part = descr.split(" - ")[0]
            end_orf = end_orf_part.split("[")[1]
        else:
            start_orf_part = descr.split(" - ")[0]
            start_orf = start_orf_part.split("[")[1]
            end_orf_part = descr.split(" - ")[1]
            end_orf = end_orf_part.split("]")[0]

        end_extended_orf = int(end_orf) + 25000 # when a range passed to blastdbcmd exceeds the length of the scaffold, blastdbcmd will by default extract the whole scaffold! :-)
        if int(start_orf) <= 25000:
            start_extended_orf = 1
        else:
            start_extended_orf = int(start_orf) - 25000

        for filename in os.listdir(blastables_outdir):
            if prot_database not in filename:
                base_name = filename.split('.')[0]
                genome_db = os.path.join(blastables_outdir, base_name)
                outscaffold = f"{genome_dir}/{res_scaf_orf}.fasta"

                success, scaffold_length = run_blastdbcmd(id, genome_db, outscaffold, f"{start_extended_orf}-{end_extended_orf}")
                if success:
                    break

print('Scaffolds of interest:', scaffolds_of_interest)


# In[71]:


# ### list containing IDs of scaffolds (not orfs) - input for blastn


#set protoze potrebujeme scaffold pouzit jako seqidlist v ramci blastn jen jednou
set_scaffolds_of_interest = set(scaffolds_of_interest)
len(set_scaffolds_of_interest)


# ted potrebujeme srovnat scaffolds_of_interest s genome_scaffold_ids a pak to pouzijeme pro seqidlist
# ziskame dict: genomes a jejich filtered scaffolds (vice virovych orfu na jednom skafoldu -> redundance)
real_genome_scaffold_ids = {}

for genome, ids in genome_scaffold_orfs_descr.items():
    truncated_ids = [id.split(';')[0] for id in ids] 
    real_genome_scaffold_ids[genome] = truncated_ids

filtered_genome_scaffold_ids = {}

for genome, ids in real_genome_scaffold_ids.items():
    for scaffold_id in ids:
        truncated_id = scaffold_id.split('_')[:-1] 
        truncated_id = '_'.join(truncated_id)
        if truncated_id in set_scaffolds_of_interest:
            if genome not in filtered_genome_scaffold_ids:
                filtered_genome_scaffold_ids[genome] = []  # initialize the list if it doesn't exist
            filtered_genome_scaffold_ids[genome].append(scaffold_id)  

print(filtered_genome_scaffold_ids)



sum(len(id) for id in filtered_genome_scaffold_ids.values())


# ## BLASTn to identify viral repeats 

def blastn_repeats(query, output, database, seqid, outfmt):
    command = ["blastn",
               "-query", str(query), 
               "-out", str(output),
               "-db", str(database), 
               "-strand", "minus",
              # "-dust", "no",
            #   "-word_size", "9",
               "-seqidlist", str(seqid), 
               "-num_threads", "5", # zmenit na clusteru
               "-outfmt", str(outfmt)]
    subprocess.run(command, check=True) 


output_dir = "/home/majnusova/all/projects/plv/paratrimastix/blastn_tirs/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
blastables_outdir = "/home/majnusova/all/projects/plv/paratrimastix/blastables/"
if not os.path.exists(blastables_outdir):
    os.makedirs(blastables_outdir)
path_to_input = "/home/majnusova/all/projects/plv/paratrimastix/eve_outscaffolds/"
if not os.path.exists(path_to_input):
    os.makedirs(path_to_input)

for genome, scaffold_ids in filtered_genome_scaffold_ids.items():
    database = f"{blastables_outdir}/{genome}" 
    outdir = os.path.join(output_dir, genome)
    os.makedirs(outdir, exist_ok=True)
    
    for scaffold_id in scaffold_ids:  
        query_file = f"{scaffold_id}.fasta" # scaffold_id = NODE_135_length_281329_cov_13.280030_127 -> NODE_135_length_281329_cov_13.280030_127.fasta
        query = os.path.join(path_to_input, genome, query_file) #/home/majnusova/all/projects/plv/data/eve_outscaffolds/NODE_135_length_281329_cov_13.280030_127.fasta
        seqidlist_file = '_'.join(scaffold_id.split('_')[:-1]) #NODE_135_length_281329_cov_13.280030
        seqidlist_path = os.path.join(outdir, f"{seqidlist_file}_seqid.txt") #/home/majnusova/all/projects/plv/data/alias_test/NODE_135_length_281329_cov_13.280030
        
        with open(seqidlist_path, 'w') as file:
            file.write(seqidlist_file + '\n')  

        output_file = f"{scaffold_id}_repeats.txt"
        output = os.path.join(outdir, output_file)
        output_file2 = f"{scaffold_id}_repeats_human.txt"
        output2 = os.path.join(outdir, output_file2)

    # Checking if the query file exists to avoid the error!
        if os.path.exists(query):
            # outfmt 6 output
            blastn_repeats(query, output, database, seqidlist_path, 5) # should be 5
            # human-readable output
            blastn_repeats(query, output2, database, seqidlist_path, 0)
        else:
            print(f"Query file {query} does not exist, skipping this scaffold ID.")


# ## Searching for suitable repeats in the blast outputs

hsp_info = defaultdict(list)
hsp_tirs_temp = defaultdict(list)  
hsp_tirs_filtered = defaultdict(list) 

blastout_dir = "/home/majnusova/all/projects/plv/paratrimastix/blastn_tirs/"
if not os.path.exists(blastout_dir):
    os.makedirs(blastout_dir)
    
for genome in filtered_genome_scaffold_ids.keys():
    fullpath = os.path.join(blastout_dir, genome)
    for blastout in os.listdir(fullpath):
        blastout_path = os.path.join(fullpath, blastout)
        if blastout.endswith("repeats.txt"):
            with open(blastout_path) as blast_handle:
                blast_records = NCBIXML.parse(blast_handle)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            if 200 <= len(hsp.match) <= 1200:
                                key = (genome, hsp.bits, hsp.score, hsp.expect, len(hsp.match))  
                                hsp_info[key].append(hsp)
                                tir = hsp.query.replace("-", "")
                                hsp_tirs_temp[(key, blastout)].append((hsp.query_start, tir))  

for (key, blastout), tirs_with_positions in hsp_tirs_temp.items():
    if len(tirs_with_positions) == 2: #interested in only two tirs - more repeats are usually connected to the problems with repeats assembling etc.
        tirs_with_positions.sort() # ascending sorting to avoid negative results
        if tirs_with_positions[1][0] - tirs_with_positions[0][0] >= 17000: # tuple(position, tir), druha[1]-prvni[0] tuple v listu, [1][0] = druha tuple[1] a jeji query_start[0]
            hsp_tirs_filtered[(key, blastout)] = tirs_with_positions


for (key, blastout), tirs in hsp_tirs_filtered.items():
    genome = key[0]
    print(f"genome: {genome}, node: {blastout}, params: {key}, number of TIRs: {len(tirs)}")
    for position, tir in tirs:
        print(f"    start_position: {position}, TIR: {tir}")


# Create a dictionary to store the count of blastout occurrences
blastout_count = {}

# Count the occurrences of each blastout value
for (key, blastout), _ in hsp_tirs_filtered.items():
    if blastout in blastout_count:
        blastout_count[blastout] += 1
    else:
        blastout_count[blastout] = 1

# Filter out non-unique blastout entries from hsp_tirs_filtered - when scaffolds contains multiple pairs of repeats, it's shit 
hsp_tirs_filtered = {k: v for k, v in hsp_tirs_filtered.items() if blastout_count[k[1]] == 1}

# Print the filtered dictionary
for (key, blastout), tirs in hsp_tirs_filtered.items():
    genome = key[0]
    print(f"genome: {genome}, node: {blastout}, params: {key}, number of TIRs: {len(tirs)}")
    for position, tir in tirs:
        print(f"    start_position: {position}, TIR: {tir}")


# ## Turning extracted scaffolds into single line ones to make it possible to search for TIRs in them

with open("/home/majnusova/all/projects/plv/data/filt_teststartpositions.txt", "w") as f:
    for (key, blastout), tirs in hsp_tirs_filtered.items():
        f.write(f"{key},{blastout}\n{tirs}\n")


scaffolds_dir = "/home/majnusova/all/projects/plv/paratrimastix/eve_outscaffolds/"
singleline_dir = "/home/majnusova/all/projects/plv/paratrimastix/eve_singleline_scaffolds/"
if not os.path.exists(singleline_dir):
    os.makedirs(singleline_dir)

for genome, ids in genome_scaffold_orfs_descr.items():
    absolute_path = os.path.join(scaffolds_dir, genome)
    output_genome_dir = os.path.join(singleline_dir, genome)  
    if not os.path.exists(output_genome_dir):
        os.makedirs(output_genome_dir)
    
    for outscaffold in os.listdir(absolute_path):
        outscaffold_path = os.path.join(absolute_path, outscaffold) 
        with open(outscaffold_path, "r") as infile:
            lines = infile.readlines()
        header = lines[0]  # Fasta header
        joined_lines = ''.join(line.strip() for line in lines[1:])  # joining all lines into one except for the first one (header) 
        
        output_file_path = os.path.join(output_genome_dir, outscaffold)
        
        with open(output_file_path, "w") as outfile:
            outfile.write(header) 
            outfile.write(joined_lines)


# ## Searching for TIRs (saved in hsp_tirs) in single-line scaffolds



singleline_dir = "/home/majnusova/all/projects/plv/paratrimastix/eve_singleline_scaffolds/"
tir_dir = "/home/majnusova/all/projects/plv/paratrimastix/eves_tirs/"

if not os.path.exists(tir_dir):
    os.makedirs(tir_dir)

for key_blastout, repeats in hsp_tirs_filtered.items():
    key, blastout = key_blastout[0], key_blastout[1]
    genome = key[0]

    genome_dir = os.path.join(tir_dir, genome)
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)
    fasta_filename = blastout.replace("_repeats.txt", ".fasta")
    absolute_path = os.path.join(singleline_dir, genome, fasta_filename)

    if os.path.exists(absolute_path):
        with open(absolute_path, "r") as infile:
            scaffold = infile.read()
            modified = False

            for position, repeat_sequence in repeats:
                if repeat_sequence in scaffold:
                    scaffold = scaffold.replace(repeat_sequence, f"\n\n{repeat_sequence}\n\n")
                    modified = True

            if modified:
                lines = scaffold.split("\n")
                new_scaffold = lines[0]  
                for line in lines[1:]:
                    if line.strip():  
                        new_line = "\n".join([line[i:i+80] for i in range(0, len(line), 80)])  
                        new_scaffold += f"\n{new_line}\n"
                    else:
                        new_scaffold += "\n"

                output_path = os.path.join(genome_dir, fasta_filename)
                with open(output_path, "w") as outfile:
                    outfile.write(new_scaffold)
    else:
        print(f"File not found: {absolute_path}")


# ## Extracting sequence btw TIRs (gene content) and extracting ORFs using getorf + appending it to the original file
tir_dir = "/home/majnusova/all/projects/plv/paratrimastix/eves_tirs/"
orfs_dir = "/home/majnusova/all/projects/plv/paratrimastix/orfs/"

if not os.path.exists(orfs_dir):
    os.makedirs(orfs_dir)

def run_getorf(inf, outf):
    command = ["getorf", "-sequence", inf, "-outseq", outf, "-table", "0", "-minsize", "450", "-find", "1"]
    subprocess.run(command, check=True)

for folder in os.listdir(tir_dir):
    folder_path = os.path.join(tir_dir, folder)
    output_folder_path = os.path.join(orfs_dir, folder)

    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)

    for file in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file)
        with open(file_path, "r") as infile:
            scaffold = infile.read()
            sections = scaffold.split("\n\n\n")
            if len(sections) > 2:
                sequence_between_tirs = sections[2].replace("\n", "")
                orf_file_name = os.path.basename(file_path)
                orf_file_path = os.path.join(output_folder_path, f"nucl_{orf_file_name}")

                with open(orf_file_path, "w") as fileorf:
                    fileorf.write(f">{orf_file_name}\n{sequence_between_tirs}\n")

                translated_orf_file_path = os.path.join(output_folder_path, orf_file_name)
                run_getorf(orf_file_path, translated_orf_file_path)

                with open(file_path, "a") as original_file:
                    with open(translated_orf_file_path, "r") as translated_file:
                        original_file.write(translated_file.read())


