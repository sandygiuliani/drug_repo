# config.py - settings for drug_repo.py


############################################################################
### PERSONAL_INFO
############################################################################
# your name
your_name = "Sandra"

# email (for NCBI Expasy)
your_email = "sandraxgiuliani@gmail.com"
############################################################################




############################################################################
### TAXONOMY
############################################################################
# define taxa as the list of taxonomy identifiers we are interested in
# e.g. SCHMA (S. Mansoni), SCHHA (S. haematobium), SCHJA (S. japonicum)
taxa = ['SCHMA']

species_lst = []

species = ''

taxa_to_species = {'SCHMA': 'S. mansoni',
					'SCHHA':'S. haematobium','SCHJA':'S. japonicum'}

# below here make string with names of species
for item in taxa:
	species_lst.append(taxa_to_species[item])


if len(species_lst) == 1:
	species = str(species_lst[0])

if len(species_lst) == 2:
	species = (species_lst[0] + ' and ' + species_lst[1])

if len(species_lst) > 2:
	for i in range(0,(len(species_lst)-2)):
		species = (species + species_lst[i] + ', ') 

	species = (species + species_lst[len(species_lst)-2] + ' and ')
	species = (species + species_lst[len(species_lst)-1])

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

############################################################################




############################################################################
### CHEMBL FILTERING SETTINGS
############################################################################
# define list of clinical phases we are interested
# eg. '4', '3', '' (empty string for the unknown phase)
chembl_phases = ['4']

# define molecule types we are interested in
chembl_mol_type = ['Synthetic Small Molecule']
############################################################################




############################################################################
### CLUSTERING SETTINGS
############################################################################
# similarity threshold for clustering
sim_threshold = 0.9
# check if float
#logger.info(isinstance(SIM_THRESHOLD, float))
############################################################################




############################################################################
### REPOSITIONING CANDIDATE
############################################################################
# repositioning candidate you wish to examine
# put CHEMBL of DB ID eg 'CHEMBL98'
repo_candidate = 'CHEMBL941'


############################################################################



############################################################################
### REGEX
############################################################################
# import regex
import re

# format of CATH domain eg '4.10.400.10'
cath_format = re.compile('.*\..*\..*\..*')

# regular expression for string containing at least one dash
contains_dash = re.compile('.*-.*')

# regular expression for string containing at least one '#'
contains_comment = re.compile('#.*')

# chembl id format
chembl_format = re.compile('CHEMBL.*')

# drugbank id format
drugbank_format = re.compile('DB.*')

# starts with colon
starts_colon = re.compile(':.*')
############################################################################
