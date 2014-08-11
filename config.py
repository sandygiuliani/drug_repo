# Copyright 2014 Sandra Giuliani
# config.py

# Configuration file for drug_repo.py


############################################################################
### PERSONAL_INFO
############################################################################
# what is your name?
your_name = "Sandra"

# what is your email? (for NCBI Expasy and T-coffee)
your_email = "sandraxgiuliani@gmail.com"
############################################################################




############################################################################
### PIPELINE STEPS
############################################################################
# define which steps of the pipeline you wish to run
# integer between 0 and 10
# eg steps = 6 will run all steps up to (and including) 6
steps = 8
############################################################################




############################################################################
### TAXONOMY
############################################################################
# define list of taxa ids you are interested in
# eg ['SCHMA', 'SCHHA', 'SCHJA']
taxa = ['SCHMA']

# to identify a specific species, look up species name in speclist.txt 
# to find the mnemonic code
# e.g. Schistosoma
# SCHMA (S. Mansoni), SCHHA (S. haematobium), SCHJA (S. japonicum)
# e.g Trypanosoma
# TRYB2 = Trypanosoma brucei brucei (strain 927/4 GUTat10.1)
# TRYB9 = Trypanosoma brucei gambiense (strain MHOM/CI/86/DAL972)
# TRYBB = Trypanosoma brucei brucei
# TRYBG = Trypanosoma brucei gambiense
# TRYBR = Trypanosoma brucei rhodesiense
# TRYCC = Trypanosoma cruzi (strain CL Brener)
# TRYCI = Trypanosoma congolense (strain IL3000)
# TRYCO = Trypanosoma congolense
# TRYCR = Trypanosoma cruzi
# TRYEQ = Trypanosoma equiperdum
# TRYEV = Trypanosoma evansi
# e.g. plasmodium (there are many others!)
# PLAF1 E   57265: N=Plasmodium falciparum (isolate 311)
############################################################################




############################################################################
### PATHS
############################################################################
# path to archindex binary
# old path "./../archSchema/bin/archindex" still valid on mac
# new path on linux machine "./../Arch/archindex"
archindex_path = "./../Arch/archindex"

# absolute path to SMSD directory (where SMSD.sh is)
# 1.5.1 - first version I have used (from sourceforge)
# 1.6 - version sent by Asad that should handle multiple sdf and keep ids
smsd_path = "/home/sandra/SMSD1.6"
############################################################################




############################################################################
### INPUT_FILES
############################################################################
# drug file from ChEMBL ('Browse drugs') 'chembl_drugs.txt'
# number of drugs should be 10406
chembl_input = 'chembl_drugs.txt'

# define CHEMBL_TARGETS as the target file from ChEMBL ('Browse drug targets')
# number of drugs associated with targets should be 2007
chembl_targets = 'chembl_drugtargets.txt'

# define CHEMBL_UNIPROT as the chemblID/uniprot mapping file
chembl_uniprot = 'chembl_uniprot_mapping.txt'

# define DRUGBANK_INPUT as the DrugBank Drug Target Identifiers
# either: all_target_ids_all.csv (all drugs, 4,026 entries),
# or: small_molecule_target_ids_all.csv (small molecule drugs, 3,899 entries)
drugbank_input = 'small_molecule_target_ids_all.csv'

# define sdf file with drugbank drugs (contains smiles)
drugbank_sdf = 'all.sdf'

# uniprot to pdb csv mapping file
# if necessary, uniprot_pdb.tsv (tsv version) can be retrieved
uniprot_pdb = "uniprot_pdb.csv"

# pdb to lig mapping file
pdb_lig = "lig_pairs.lst"

# pointless het groups
pointless_het = "pointless_het.csv"

# chemical component smiles dictionary
cc_smi = "Components-smiles-oe.smi"

# location of the species codes to species names mapping file 
spec_list = 'speclist.txt'
############################################################################




############################################################################
### OUTPUT_FILES
############################################################################
# define names of output files, they will be overwritten every time
# if you do not want that to happen, add a timestamp to the file names
# 'dr' stands for drug repositioning

# log
log_name = 'dr_log.log'



############################################################################




############################################################################
### CHEMBL FILTERING SETTINGS
############################################################################
# define list of clinical phases you are interested in 
# (only applies to ChEMBL set)
# eg. '4', '3', '' (empty string for the unknown phase)
chembl_phases = ['4']

# define molecule types you are interested in
chembl_mol_type = ['Synthetic Small Molecule']
############################################################################




############################################################################
### CLUSTERING SETTINGS
############################################################################
# define similarity threshold for clustering
# e.g. 0.9
sim_threshold = 0.9
# check if float
#logger.info(isinstance(SIM_THRESHOLD, float))
############################################################################




############################################################################
### REPOSITIONING CANDIDATE
############################################################################
# repositioning candidate to be examined

# put CHEMBL or DB ID eg 'CHEMBL98'
repo_candidate = 'CHEMBL941'

# target number, for selecting which drug target to align to the potential
# parasite targets. 
# 0 is the first one (could be the only one), 1 the second one...
repo_target_no = 0
############################################################################




############################################################################
### HOMOLOGY MODEL
############################################################################
# number of homology models to make
model_no = 10
# alignment file - has to be in PIR format
model_align = '1d3h_schma.ali'
# template name - PDB ID of the crystal structure
model_xray = '1d3h'
# sequence to model name - arbitrary name, but has to match in the .ali file
model_seq = 'schma'
############################################################################