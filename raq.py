# # A script that finds homologous protein sequences by performing hmmsearch and control blast.

# In[1]:


from glob import glob
import re
import subprocess
from Bio import SearchIO
from Bio import SeqIO
from Bio.Blast import NCBIXML


# ### Running hmmer: HMM against the concatenated Eukprot file. The results are saved as hmmout_eukprot.txt

# In[2]:


def run_hmmsearch(database, hmm_model, output_file):
    # Build the command
    command = ["hmmsearch", "-o", output_file, hmm_model, database]
    # Execute the command
    subprocess.run(command, check=True)
# variables should be defined outside of the function to be reusable
hmm_model = "/home/users/mvladka/projects/jotnarlogs/data/raq.hmm"
database = "/home/users/mvladka/eukprot/proteins/eukprot_concatenated.fasta"
output_file = "/home/users/mvladka/projects/jotnarlogs/intermediate_files/raq_hmmout_eukprot.txt"


# In[3]:


run_hmmsearch(database, hmm_model, output_file) 


# ### Parsing hmmer outputs and filtering only those results that have e-value better than e-05 and belong to the organisms we want to map.

# In[6]:


best_hmmer_hits = [] 

hmmout_parsed = SearchIO.parse(output_file, 'hmmer3-text')
for result in hmmout_parsed: # iterating through the results of hmmsearch
    for hit in result: 
        evalue = hit.evalue  # Access evalue from the HIT object (hit vs. hsp, hit may be composed of multiple HSPs)
        if "e" in str(evalue) or float(evalue) == 0: # we wanted to set the treshold to e-05 (-05 is problematic for python -> usage of "e")
            best_hmmer_hits.append(hit.id) # appending IDs of hits with evalue better than e-05


# In[ ]:





# In[3]:


len(best_hmmer_hits)


# In[28]:


orgs_to_map = ["Homo_sapiens",
"Takifugu_rubripes",
"Ciona_intestinalis",
"Branchiostoma_floridae",
"Strongylocentrotus_purpuratus",
"Saccoglossus_kowalevskii",
"Drosophila_melanogaster",
"Daphnia_pulex",
"Caenorhabditis_elegans",
"Adineta_vaga",
"Lottia_gigantea",
"Helobdella_robusta",
"Capitella_teleta",
"Schistosoma_mansoni",
"SMU",
"Nematostella_vectensis",
"Hydra_magnipapillata",
"KII",
"Trichoplax_adhaerens",
"Mnemiopsis_leidyi",
"Amphimedon_queenslandica",
"Oscarella_carmela",
"Monosiga_brevicollis",
"Salpingoeca_rosetta",
"Capsaspora_owczarzaki",
"Pigoraptor_chileana",
"Sphaeroforma_arctica",
"Creolimax_fragrantissima",
"Amoebidium_parasiticum",
"Chromosphaera_perkinsii",
"Syssomonas_multiformis",
"Corallochytrium_limacisporum",
"Tunicaraptor_unikontum",
"Coprinopsis_cinerea",
"Cryptococcus_neoformans",
"Ustilago_maydis",
"Schizosaccharomyces_pombe",
"Saccharomyces_cerevisiae",
"Neurospora_crassa",
"Aspergillus_fumigatus",
"Phycomyces_blakesleeanus",
"Coemansia_reversa",
"Conidiobolus_coronatus",
"Olpidium_bornovanus",
"Basidiobolus_meristosporus",
"Catenaria_anguillulae",
"Amoeboradix_gromovi",
"Sanchytrium_tribonematis",
"Batrachochytrium_dendrobatidis",
"Spizellomyces_punctatus",
"Gonapodya_prolifera",
"Hyaloraphidium_curvatum",
"Paraphelidium_tribonemae",
"Rozella_allomycis",
"Paramicrosporidium_saccamoebae",
"Encephalitozoon_cuniculi",
"Fonticula_alba",
"Thecamonas_trahens",
"Amastigomonas_sp",
"Pygsuia_biforma",
"Dictyostelium_discoideum",
"Physarum_polycephalum",
"Entamoeba_histolytica",
"Mastigamoeba_balamuthi",
"Idionectes_vortex",
"Armaparvus_languidus",
"Acanthamoeba_castellanii",
"Rigifila_ramosa",
"Diphylleia_rotans",
"Mantamonas_plastica",
"Mantamonas_vickermani",
"Mantamonas_sphyraenae",
"Planomonas_micra",
"Ancyromonas_sigmoides",
"Malawimonas_sp_californiana",
"Gefionella_okellyi",
"Trichomonas_vaginalis",
"Giardia_intestinalis",
"Spironucleus_vortens",
"Spironucleus_salmonicida",
"Retortamonas_caulleryi",
"Dysnectes_brevis",
"Kipferlia_bialata",
"Chilomastix_cuspidata",
"Chilomastix_caulleryi",
"Aduncisulcus_paluster",
"Ergobibamus_cyprinoides",
"Euthynema_mutabile",
"Iotanema_spirale",
"Carpediemonas_membranifera",
"Barthelona_sp_PAP020",
"Monocercomonoides_exilis",
"Paratrimastix_pyriformis",
"Anaeramoeba_ignava",
"Anaeramoeba_flamelloides",
"Trypanosoma_brucei",
"Trypanosoma_cruzi",
"Leishmania_major",
"Bodo_saltans",
"Perkinsela_sp",
"Apiculatamorpha_spiralis",
"Euglena_gracilis",
"Naegleria_gruberi",
"Neovahlkampfia_damariscottae",
"Andalucia_godoyi",
"Agogonia_voluta",
"Ophirina_chinija",
"Arabidopsis_thaliana",
"Oryza_sativa",
"Ginkgo_biloba",
"Selaginella_moellendorffii",
"Physcomitrella_patens",
"Penium_margaritaceum",
"Klebsormidium_nitens",
"Ostreococcus_lucimarinus",
"Micromonas_commoda",
"Coccomyxa_subellipsoidea",
"Chlamydomonas_reinhardtii",
"Volvox_carteri",
"Cyanophora_paradoxa",
"Gloeochaete_wittrockiana",
"Cyanidioschyzon_merolae",
"Galdieria_sulphuraria",
"Porphyridium_purpureum",
"Chondrus_crispus",
"Rhodelphis_limneticus",
"Rhodelphis_marinus",
"Picomonas_judraskeda",
"Picobiliphyte_sp._MS584-11",
"Picobiliphyte_sp._MS584-22",
"Picobiliphyte_sp._MS584-5",
"Picobiliphyte_sp._MS609-66",
"Guillardia_theta",
"Goniomonas_avonlea",
"Hemiarma_marina",
"Palpitomonas_bilix",
"Microheliella_maris",
"Emiliania_huxleyi",
"Pavlova_pinguis",
"Raineriophrys_erinaceoides",
"Ancoracysta_twista",
"Colponemidia_sp_Colp-10",
"Tetrahymena_thermophila",
"Paramecium_tetraurelia",
"Chromera_velia",
"Cryptosporidium_parvum",
"Toxoplasma_gondii",
"Plasmodium_falciparum",
"Theileria_annulata",
"Perkinsus_marinus",
"Symbiodinium_minutum",
"Aurantiochytrium_limacinum",
"Aplanochytrium_kerguelense",
"Blastocystis_sp_subtype7",
"Cafeteria_burkhardae",
"Phytophthora_sojae",
"Hyaloperonospora_arabidopsidis",
"Thalassiosira_pseudonana",
"Phaeodactylum_tricornutum",
"Aureococcus_anophagefferens",
"Microchloropsis_gaditana",
"Vischeria_C74",
"Ectocarpus_siliculosus",
"Bigelowiella_natans",
"Paulinella_micropora",
"Plasmodiophora_brassicae",
"Reticulomyxa_filosa",
"Telonema_subtile",
"Hemimastix_kukwesjijk",
"Spironema_sp_BW2"]


# In[40]:


wanted_sequences_orgs = set()
str_best_hmmer_hits = str(best_hmmer_hits)

for org in orgs_to_map:
    if org.strip() in str_best_hmmer_hits: # faster than str(best_hmmer_hits)
        for selected_org in best_hmmer_hits: 
            if org.strip() in selected_org:
                wanted_sequences_orgs.add(selected_org)

wanted_sequences_orgs


# In[42]:


len(wanted_sequences_orgs)


# In[33]:


# print(dir(hit), end="")


# #### Saving the IDs of hits (e-05, selected orgs) into a file "hmmer_05_hits.txt"

# In[39]:


with open('/home/users/mvladka/projects/jotnarlogs/intermediate_files/raq_hmmer_05_hits.txt', 'w') as outf:
    for hit_id in wanted_sequences_orgs:
        outf.write(f"{hit_id}\n")


# ## Extracting fasta sequences based on the HMMER results using the blastdbcmd command

# In[8]:


"""def makeblastdb(input_fas, blast_db):
    command = ['makeblastdb','-in', input_fas,'-dbtype', 'prot','-out', blast_db, '-parse_seqids']
    subprocess.run(command, check=True)

input_fas = '/home/users/mvladka/eukprot/proteins/eukprot_concatenated.fasta' 
blast_db = '/home/users/mvladka/eukprot/proteins/eukprot_concatenated' #uz je hotovo

makeblastdb(input_fas, blast_db)"""


# In[10]:


def extract_sequences(database, entry_batch_file, output_fasta):
    extract_sequences_command = [
        'blastdbcmd',
        '-db', database,
        '-entry_batch', entry_batch_file,
        '-outfmt', '%f',
        '-out', output_fasta
    ]
    subprocess.run(extract_sequences_command, check=True)
database = '/home/users/mvladka/eukprot/proteins/eukprot_concatenated'
entry_batch_file1 = '/home/users/mvladka/projects/jotnarlogs/intermediate_files/raq_hmmer_05_hits.txt'
output_fasta = '/home/users/mvladka/projects/jotnarlogs/intermediate_files/raq_candidates.fas'


# In[11]:


extract_sequences(database, entry_batch_file1, output_fasta)


# ### Running control BLASTP: candidate seqs x RABs database.

# In[12]:


blastdb = '/home/users/mvladka/projects/jotnarlogs/data/rabs_final.txt'
query_file = '/home/users/mvladka/projects/jotnarlogs/intermediate_files/raq_candidates.fas'
blastout_xml = '/home/users/mvladka/projects/jotnarlogs/intermediate_files/raq_control_blast_result.xml'
blastout_text = '/home/users/mvladka/projects/jotnarlogs/intermediate_files/raq_control_blast_result.txt'

# Run BLAST with XML output
cmd_xml = ['blastp', '-query', query_file, '-db', blastdb, '-outfmt', '5', '-out', blastout_xml]
subprocess.run(cmd_xml, check=True, text=True)

# Run BLAST with human-readable output 
cmd_text = ['blastp', '-query', query_file, '-db', blastdb, '-out', blastout_text]
subprocess.run(cmd_text, check=True, text=True)


# ### Checking the results of control blast and looking for hits containing RAQ/RAW...

# In[15]:


blast_results = '/home/users/mvladka/projects/jotnarlogs/intermediate_files/raq_control_blast_result.xml'
wanted_sequences = set()

blast_records = NCBIXML.parse(open(blast_results))
for blast_record in blast_records:
    query_sequence = blast_record.query # .query = 'EP00134_Neurospora_crassa_P007575 XP_961711.1 GTP-binding protein GTR1 [Neurospora crassa OR74A]'

    alignment_count = 0  #counter for the number of alignments processed (nekdy je raw az uplne dole a jen jednou - nebude to ono)

    for alignment in blast_record.alignments:
        alignment_count += 1
        if alignment_count > 20:  # Check if more than 20 alignments have been processed
            break

        hit_descr = alignment.title
        if "raq" in hit_descr.lower() or "raxx1" in hit_descr.lower() or "raxx4" in hit_descr.lower(): #prepsat pak na RAQ apod. 
            seqid = query_sequence.split(" ")[0].strip()
            wanted_sequences.add(seqid)

# In[17]:


# Save the unique hit IDs to a file
with open('/home/users/mvladka/projects/jotnarlogs/intermediate_files/raq_control_blast_ids.txt', 'w') as outf:
    for raq_id in wanted_sequences:
        outf.write(f"{raq_id}\n")


# ### Extracting the putative RAQ, RAW seqs into a file that should be checked manually.

# In[19]:


database = '/home/users/mvladka/eukprot/proteins/eukprot_concatenated' 
entry_batch_file2 = '/home/users/mvladka/projects/jotnarlogs/intermediate_files/raq_control_blast_ids.txt'
output_fasta2 = '/home/users/mvladka/projects/jotnarlogs/results/raq_final_candidates.fas'

extract_sequences(database, entry_batch_file2, output_fasta2)
