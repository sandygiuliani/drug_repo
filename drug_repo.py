# Copyright 2014 Sandra Giuliani
# drug_repo.py

# Mapping of known drugs from ChEMBL and DrugBank to their targets/domain
# architecture to identify suitable drug repositioning candidates for
# schistosomiasis.
# Please see README.md for more info.




############################################################################
### import modules, set up logger
############################################################################

# import os (old system)
import os
# os system does not work very well, it is deprecated, use subprocess
#os.system('./../archSchema/bin/archindex -u B6DTB2 > test2.txt')

# import os.path for checking if files exist
import os.path

# import pickle
import cPickle as pickle

# import csv for comma-sep files
import csv

# import subprocess for executing command line
import subprocess

# import itertools for flatten out lists
import itertools

# import other modules
import sys, re, string, fnmatch, shutil

# import Biopython Entrez
from Bio import Entrez
# tell NCBI who I am
Entrez.email = "sandraxgiuliani@gmail.com"

# import SeqIO
from Bio import SeqIO

# list available databases
#handle = Entrez.einfo()
#logger.info(handle.read())

# import expasy for access protein sequences
from Bio import ExPASy

# import swissprot for parsing swissprot plain text files
from Bio import SwissProt

# set up log
import logging
# set up log file to write to, it will be overwritten every time ('w' mode)
# leave this level setting to DEBUG
logging.basicConfig(filename='log_drug_repo.log', filemode='w',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)
# leave this level setting to DEBUG
logger.setLevel(logging.DEBUG)
# create console handler
ch = logging.StreamHandler()
# CHANGE THIS TO TUNE LOGGING LEVEL from DEBUG/INFO/WARNING
ch.setLevel(logging.DEBUG)
# create formatter, you can add:
# '%(levelname)s' for level eg DEBUG, INFO..
# '%(name)s' for level name, eg __main__ in the log
formatter = logging.Formatter('%(asctime)s - %(message)s')
# add formatter to ch
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)
############################################################################




############################################################################
### define global variables
############################################################################

# define CHEMBL_INPUT as the drug file from ChEMBL ('Browse drugs')
# 'chembl_drugs.txt'
# number of drugs should be 10406
CHEMBL_INPUT = 'chembl_drugs.txt'

# define CHEMBL_TARGETS as the target file from ChEMBL ('Browse drug targets')
# number of drugs associated with targets should be 2007
CHEMBL_TARGETS = 'chembl_drugtargets.txt'

# define CHEMBL_UNIPROT as the chemblID/uniprot mapping file
CHEMBL_UNIPROT = 'chembl_uniprot_mapping.txt'

# define list of clinical phases we are interested
# eg. '4', '3', '' (empty string for the unknown phase)
CHEMBL_CLINICAL_PHASES = ['4']

# define molecule types we are interested in
CHEMBL_MOL_TYPE = ['Synthetic Small Molecule']

# define DRUGBANK_INPUT as the DrugBank Drug Target Identifiers
# either: all_target_ids_all.csv (all drugs, 4,026 entries),
# or: small_molecule_target_ids_all.csv (small molecule drugs, 3,899 entries)
DRUGBANK_INPUT = 'small_molecule_target_ids_all.csv'

# define TAXA as the list of taxonomy identifiers we are interested in
# e.g. SCHMA (S. Mansoni), SCHHA (S. haematobium), SCHJA (S. japonicum)
TAXA = ['SCHMA', 'SCHHA', 'SCHJA']

# format of CATH domain eg '4.10.400.10'
CATH_FORMAT = re.compile('.*\..*\..*\..*')

# path to archindex binary
# old path "./../archSchema/bin/archindex" still valid on mac
# new path on linux machine "./../Arch/archindex"
ARCHINDEX_PATH = "./../Arch/archindex"
############################################################################




############################################################################
### HEADER_COUNT HELPER FUNCTION
############################################################################
# find specific header in line of comma or tab (or other) separated
# file and return column number of the header
def header_count(line, separator, header):
  '''
  Read a string (line) separated by separator
  and return column number of a specific string (header)
  '''
  #set counter to 0
  col_number = 0
  # loop until word matches the header
  while not line.split(separator)[col_number] == header:
    col_number = col_number + 1
  # return column number
  return col_number
############################################################################



############################################################################
### FILE_TO_LINES HELPER FUNCTION
############################################################################
# read file using readlines approach and return the lines
def file_to_lines(text_file):
  try:
    input_file = open(text_file, 'r')
    lines = input_file.readlines()
    input_file.close()
    #logger.debug('the try/expect loop is working')
    return lines
  except IOError:
    logger.error('The file ' + text_file + ' cannot be found' +
                 ' in the current directory!')
    #this logger exception will print the whole thing
    #logger.exception('whaaaaat we have error')
    logger.warning('The program is aborted.')
    # exit python
    sys.exit()
############################################################################




############################################################################
### SWAP_DIC HELPER FUNCTION
############################################################################
# read tab-separated mapping file with header and return dictionary with
# second column as key and first column as values - created for the
# chemblID uniprot mapping file
def swap_dic(tab_file):
  lines = file_to_lines(tab_file)
  swap_dictionary = {}
  counter = 0
  # iterate over lines
  for i in range(1,len(lines)):
    # split tab
    splitline = lines[i].split("\t")
    counter = counter + 1
    #logger.debug(splitline[0])
    # create dictionary, stripping the carriage return
    swap_dictionary[splitline[1].rstrip('\r\n')] = (splitline[0])
  #logger.debug(swap_dictionary)
  return swap_dictionary
############################################################################



############################################################################
### PROCESS_CHEMBL FUNCTION
############################################################################

# in the end we need a dictionary of chembl ids vs uniprot ids
def process_chembl():
  '''read chembl drug input file, filter, get various info and return
  list of uniprot ids'''
  # open chembldrugs.txt for reading
  lines = file_to_lines(CHEMBL_INPUT)
  logger.info('We are now analysing the ChEMBL input file ' +
              str(CHEMBL_INPUT) + ', containing a total of '
              + str(len(lines)-1) + ' drugs.')


  ### ANALYSE HEADERS AND OBTAIN COLUMN NUMBERS

  # get the headers
  headers = lines[0]
  # remove duplicate spaces
  #headersnospace = " ".join(headers.split())
  # print headers list in lowercase
  # print('The headers are: '+headersnospace.lower())

  col_phase = header_count(headers, "\t", "DEVELOPMENT_PHASE")
  col_type = header_count(headers, "\t", "DRUG_TYPE")
  col_chemblid = header_count(headers, "\t", "CHEMBL_ID")

  # count the drugs in the 4+1 classes of clinical phase
  # set counters for clinical phases to zero
  phase1 = 0
  phase2 = 0
  phase3 = 0
  phase4 = 0
  phase_unknown =0
  # iterate over rows, excluding the header row
  for i in range(1,len(lines)):
    # tab separate each row
    rowsplit = lines[i].split("\t")
    if (rowsplit[col_phase] == '1'):
      phase1 = phase1 + 1
    elif (rowsplit[col_phase] == '2'):
      phase2 = phase2 + 1
    elif (rowsplit[col_phase] == '3'):
      phase3 = phase3 + 1
    elif (rowsplit[col_phase] == '4'):
      phase4 = phase4 + 1
    else:
      phase_unknown = phase_unknown + 1

  logger.info('The number of drugs in phase 1 is: ' + str(phase1) +
              '; in phase 2 is: ' + str(phase2) + '; in phase 3 is: ' +
              str(phase3)+'; in phase 4 is: ' + str(phase4) +
              '; in unknown phase is: ' + str(phase_unknown) + '.')
  ###


  ### CLINICAL PHASE AND MOLECULAR TYPE FILTER
  # IN THE END WE OBTAIN CHEMBL_FILT_LIST THAT HAD THE FILTERED CHEMBL IDS
  # possibly merge THEM TOGETHER!

  # set up empty list to append lines to
  stripped = []
  # look over lines, excluding header
  for y in range(1,len(lines)):
    # tab separate each row
    rowsplit2 = lines[y].split("\t")
    # check if they are in clinical phase we are intereste in
    if (rowsplit2[col_phase] in CHEMBL_CLINICAL_PHASES):
      # append the stripped lines to the list
      stripped.append(lines[y])

  # empty list in which to store filtered chembl drug ids
  chembl_filt_list = []
  # loop over, note here there is no header
  for line in stripped:
    # tab separate
    rowsplit3 = line.split("\t")
    # check if they are of the desired drug type
    if rowsplit3[col_type] in CHEMBL_MOL_TYPE:
      #logger.debug(rowsplit3[col_chemblid])

      # make list of chembl ids we are interested in
      chembl_filt_list.append(rowsplit3[col_chemblid])

  logger.info('We have filtered the entries in clinical phases ' +
              str(CHEMBL_CLINICAL_PHASES) + ' and of molecule type ' +
              str(CHEMBL_MOL_TYPE) + ', to obtain ' +
              str(len(chembl_filt_list)) + ' drugs.')
  ###


  # CREATE DICTIONARY CHEMBL_TARGET_DRUG_DIC IN WHICH TO STORE
  # {CHEMBL TARGET IDS: (LIST OF CHEMBL DRUG IDS)}

  # open the drug targets chembl file and get lines
  drug_targ = file_to_lines(CHEMBL_TARGETS)
  # get column number for two headers we want (chembl ids for mol and targets)
  col_mol_id = header_count(drug_targ[0], '\t', 'MOLECULE_CHEMBL_ID')
  col_targ_id = header_count(drug_targ[0], '\t', 'TARGET_CHEMBL_ID')
  #logger.debug(col_targ_id)

  # empty dictionary
  chembl_target_drug_dic = {}

  # list of target chemblids that refer to the drugs we are interested in
  # ie. small molecules, late clinical stages
  for i in range(1,len(drug_targ)):
    rowsplit = drug_targ[i].split("\t")
    #logger.debug(rowsplit)

    # get the chembl drug and target id values for the row
    chembl_drug_id = rowsplit[col_mol_id]
    chembl_target_id = rowsplit[col_targ_id]

    # only proceed if the target id is not an empty field!
    if chembl_target_id  != "":

      # check if the molecule chembl id is one of the drugs we want
      if chembl_drug_id in chembl_filt_list:

        # check if the target is is already in the dictionary, and extend
        if chembl_target_id in chembl_target_drug_dic:
          # try to append stuff to dictionary
          chembl_target_drug_dic[chembl_target_id].append(chembl_drug_id)

        # else create the new entry in the dictionary
        else:
          # create empty list
          chembl_target_drug_dic[chembl_target_id] = []
          # append to empty list the value
          chembl_target_drug_dic[chembl_target_id].append(chembl_drug_id)

  # test to see if there are duplicates drugs
  # positive: tot 2007 but 1596 unique
  #logger.debug(len(test_duplo))
  #logger.debug(len(list(set(test_duplo))))

  # number of drugs we could associate with chembl target
  # this will be a list of lists, so we need to flatten it out
  drug_ids = list(itertools.chain(*list(chembl_target_drug_dic.values())))
  #logger.debug(len(drug_ids))

  # eliminate duplicates
  drug_ids = list(set(drug_ids))
  #logger.debug(drug_ids)

  ###



  # create dictionary from the chembl/uniprot mapping file
  # the dictionary will be {'chemblID1':'uniprotid1', etc..}
  # the dictionary has to be this way because more than one chembl id
  # can point to the same uniprot id
  chembl_uniprot_map_dic = swap_dic(CHEMBL_UNIPROT)
  #logger.debug(chembl_uniprot_map_dic)

  # compare length dictionary with unique uniprot id values!
  #logger.debug(len(chembl_uniprot_map_dic))
  #logger.debug(len(list(set(chembl_uniprot_map_dic.values()))))

  
  # empty dictionary,it will store {drug chembl id : [list of uniprot]}
  chembl_dic = {}

  # new method to create swap dictionary instead!
  # loop over unique list of drug ids
  for drug in drug_ids:
    # empty target list in which to store target ids
    target_list = []
    # loop over dictionary of chembl target ids and drug chembl ids
    for key in chembl_target_drug_dic:
      # check if drug is in the list
      if drug in chembl_target_drug_dic[key]:
        # append chembl id to list
        target_list.append(key)
    #logger.debug(target_list)

    # empty list for uniprot values
    uniprot_list = []
    # loop over the list of chembl ids we have just generated
    for chembl_id in target_list:
      # check if this id is in the mapping dictionary
      if chembl_id in chembl_uniprot_map_dic:
        # append uniprot value to list
        uniprot_list.append(chembl_uniprot_map_dic[chembl_id])

    # check the list is not empty
    if uniprot_list != []:
      chembl_dic[drug] = uniprot_list

  # flat list of uniprot codes from the dictionary
  uniprot_ids = list(itertools.chain(*list(chembl_dic.values())))
  #logger.debug(len(drug_ids))

  # eliminate duplicates
  uniprot_ids = list(set(uniprot_ids))


  # NB the len of target_ids and chembl_target_drug_dic is the same!
  logger.info(str(len(chembl_dic)) + ' ChEMBL drugs could be associated ' +
              'with ' + str(len(uniprot_ids)) + ' unique UniProt ID.')
  
  #logger.debug(len(list(itertools.chain(*list(chembl_dic.values())))))
  #logger.debug(chembl_dic)

  # return dictionary {uniprot1:(list of chembl ids)}
  return chembl_dic
############################################################################




############################################################################
### PROCESS_DRUGBANK FUNCTION
############################################################################
# process drugbank file, return a dictionary of drugs vs uniprot ids. Deal
# with duplicate values and append to list in the dictionary.
def process_drugbank():

  # open and read drug_bank input and count number
  lines = file_to_lines(DRUGBANK_INPUT)
  # number of entries (remember there is header)
  logger.info('We are now analysing the DrugBank input file ' +
              str(DRUGBANK_INPUT) + ', containing ' + str(len(lines)-1) +
              ' entries.')

  # headers are first line of the file, stripped of carriage return
  headers = lines[0].rstrip('\r\n')
  #logger.debug('The DrugBank headers are: ' + headers + '.')
  # find column number of uniprot and drugbank ids
  col_uniprot = header_count(headers, "," , "UniProt ID")
  col_drugbankids = header_count(headers, "," , "Drug IDs")
  #logger.debug(col_uniprot)
  #logger.debug(col_drugbankids)

  # empty dictionary in which to store drugbank info
  drugbank_dic = {}
  # loop over to len(lines), excluding first line (headers)
  #for i in range(1,3898):
    #logger.debug(lines[i])
    # do I need here skipinitialspace=True ??

  # read csv lines with the csv reader - deals with quotation marks etc..
  incsv = csv.reader(lines)
  # skip the first line with the headers
  next(incsv)
  # loop over each line
  # this list_check to double check we are dealing with duplicates
  #list_check = []
  # loop over each line
  for line in incsv:
    #list_check.append(line[col_uniprot])
    drug_string = line[col_drugbankids]
      # list of drubbank ids
    drug_split = drug_string.split(';')
    #ogger.debug(drug_split)
      # check if uniprot id is already in the dictionary
    if line[col_uniprot] in drugbank_dic:
      # append drug id value to the list in the dictionary
      # the '.extend' prevents the formation of a list of lists!
      drugbank_dic[line[col_uniprot]].extend(drug_split)
    else:
      # populate dictionary with the new entry
      drugbank_dic[line[col_uniprot]] = drug_split

  logger.info('The DrugBank drugs could be mapped to ' +
              str(len(drugbank_dic)) + ' unique UniProt IDs.')

  # confirm we are dealing with duplicates
  #logger.debug(len(list_check))
  #logger.debug(len(list(set(list_check))))

  # return dictionary {uniprot1:(list of drugbank ids)}
  return drugbank_dic
############################################################################




############################################################################
### UNIPROT_TO_ARCH FUNCTION
############################################################################
def uniprot_to_arch(uniprot_list,architecture):
  ''' run archindex, return dictionary of CATH/pfam domain architecture vs
   uniprot values
  '''
  if architecture == "cath":
    flag = "-cath"

  elif architecture == "pfam":
    flag = ""

  # dictionary of uniprot ids and list of correposponding architectures
  arch_dic = {}
  
  # loop over list of uniprot values
  for uniprot_id in uniprot_list:
    #list in which to store list of CATH domains for each entry
    architect_list = []
    # call archschema on the list
    subprocess.call(ARCHINDEX_PATH + " -u " + str(uniprot_id) +
                  " -maxa 1 -maxs 1 " + str(flag) +" > temp.txt", shell=True)
    # store lines
    lines = file_to_lines('temp.txt')
    #logger.debug(lines)
    for i in range(len(lines)):
      # find line that starts with parent
      if lines[i][0:7] == ':PARENT':
        # take the line after the ':PARENT' and split it
        line_split = lines[i+1].split("\t")
        #logger.debug('the line after is ' + str(line_split))


        ### for cath ###
        if architecture == "cath":

          # check if there are 'p's and get rid of them
          if "p" in line_split[2]:
            # replace p's with nothing
            line_nops = line_split[2].replace('p','')
          else:
            line_nops = line_split[2]

          # check if there are undescores
          if "_" in line_nops:
            undersc_split = line_nops.split("_")
            #logger.debug(undersc_split)

            for item in undersc_split:
              # check if the format is CATH one
              #cath_format = re.compile('.*\..*\..*\..*')
              if CATH_FORMAT.match(item):
                architect_list.append(item)

          # this is the case of just one entry, no undescores
          else:
            # check the format is CATH one
            if CATH_FORMAT.match(line_nops):
              architect_list.append(line_nops)


        ### for pfam ###
        elif architecture == "pfam":
          # check if there are dots
          if "." in line_split[2]:
            dot_split = line_split[2].split(".")
            #logger.debug(undersc_split)

            for item in dot_split:

              architect_list.append(item)

          # this is the case of just one entry, no dots
          else:
            architect_list.append(line_split[2])



    # eliminate duplicates within the list (this is for each entry!!)
    architect_list = list(set(architect_list))

    # populate the dictionary
    arch_dic[uniprot_id] = architect_list


  #logger.debug(cath_dic)

  # rm temp.txt in the end
  # this is the last temp file that overwrote the others
  subprocess.call("rm temp.txt", shell=True)

  #logger.info('We have found ' + str(len(architect_list)) +
  #            ' unique CATH domain architectures.')
  
  #logger.debug(architect_list)
  # return the list of unique domain architecture values
  return arch_dic
############################################################################




############################################################################
### ARCH_TO_UNIPROT FUNCTION
############################################################################
# run archindex, filter for TAXA and find uniprot ids
def arch_to_uniprot(arch_list,architecture):
  '''(list of str -> list of str)
  run archindex, return dictionary of uniprot from list of domain
  architecture values applying a filter for TAXA (organism)
  '''

  if architecture == "cath":
    flag = "-cath"
  elif architecture == "pfam":
    flag = ""

  # empty dictionary
  uniprot_dic = {}
  # loop over list of arch values
  for arch_id in arch_list:
    # empty list in which to store uniprot values
    uniprot_list = []
    # iterate over the taxa code list (schisto species)
    for taxa_code in TAXA:
      # call archschema on the list
      subprocess.call(ARCHINDEX_PATH + " -p " + str(arch_id) +
                    " -maxa 1 -maxs 1 " + str(flag) + " -s " + taxa_code +
                    " > temp.txt", shell=True)
      # store lines
      lines = file_to_lines('temp.txt')
      #logger.debug(temp.txt)

      # get the uniprot values and append them to uniprot_list
      for i in range(len(lines)):
        # find line that starts with parent
        if lines[i][0:7] == ':PARENT':
          # take the line after the ':PARENT' and split it
          line_split = lines[i+1].split("\t")

          uniprot_list.append(line_split[0])
    
    # remove duplicates from list
    uniprot_list = list(set(uniprot_list))

    # populate the dictionary
    uniprot_dic[arch_id] = uniprot_list

  # rm temp.txt in the end
  # this is the last temp file that overwrote the others
  subprocess.call("rm temp.txt", shell=True)

  #logger.info('Using ' + str(architecture) + 
  #            ' domain architectures we have found ' + 
   #           str(len(uniprot_list)) + 
   #           ' UniProt IDs of schistosoma proteins (taxonomic identifiers ' +
   #            str(TAXA) + ').')

  #logger.debug(uniprot_list)

  #return the dictionary
  return uniprot_dic
############################################################################




############################################################################
### MERGE_LISTS
############################################################################
# merge lists, remove duplicates and return merged list
def merge_lists(list1, list2):

  merged_list = list(set(list1)|set(list2))

  #pickle.dump(uniprot_schisto_list, open("uniprot_schisto_list.p", "wb"))

  #logger.debug('The merged UniProt values obtained from CATH and pfam ' +
  #  'are ' + str(len(merged_list)) + '.')

  #logger.debug(uniprot_schisto_list)
  return merged_list
############################################################################




############################################################################
### EXPASY_FILTER
############################################################################
# take list of uniprot, call expasy and applies filter to find the reviewed
# entries

  # entrez fetch does not seem to handle non-swissprot ids (TrEMBL)
  #prothandle = Entrez.efetch(
  #                  db="protein", id= 'P13566', rettype = "gb")
  #seq_record = SeqIO.read(prothandle, "gb")
  #logger.debug(seq_record)

def expasy_filter(uniprot_list):

  # empty list to sore filtered entries
  filtered_list = []

  # counter for reviewed entries
  rev_count = 0

  # counter for records with multiple accession numbers
  #mult_count = 0

  for entry in uniprot_list:

    handle = ExPASy.get_sprot_raw(entry)
    # swissprot read
    record = SwissProt.read(handle)
    # alternative to swissport read, should also work
    #record = SeqIO.read(handle, "swiss")

    if record.data_class == 'Reviewed':
      rev_count = rev_count + 1
      filtered_list.append(entry)

    #if len(record.accessions) > 1:
    #  mult_count = mult_count + 1
    #  logger.debug('The entry ' + entry + ' has multiple entries ' + 
    #                str(record.accessions))

    #print(record.seqinfo)
    #for ref in record.references:
    #  print(ref.authors)

  logger.info('The reviewed entries are ' + str(rev_count) + '.')

  return filtered_list
  #logger.info('The entries with multiple ids are ' + str(mult_count) + '.')

############################################################################




############################################################################
### FLATTEN_VALUES
############################################################################
# list from dictionary values and flatten it (from list of lists to simple
# list, also eliminate duplicates

def flatten_values(dic):
  # generate list and flatten it
  flat_list = list(itertools.chain(*list(dic.values())))

  # rm duplicates
  flat_list = list(set(flat_list))

  return flat_list
############################################################################




############################################################################
### CHEMBL_REPO_MAP
############################################################################
# make big map for chembl ids, for both cath and pfam arch
# format {drug1:{drug_target1:{arch1:[list os schisto uniprot], arch2:[..]},
#                 drug_target2:{..}, drug2:{........}}}

# takes chembl to uniprot dic, uniprot to cath dic, cath to schisto dic,
# uniprot to pfam dic, pfam to schisto dic

def chembl_repo_map(chembl_dic, cath_dic, schisto_cath_dic, 
                    pfam_dic, schisto_pfam_dic):
  # empty dictionary
  chembl_repo_map = {}

  # add stuff

  return chembl_repo_map
############################################################################




############################################################################
### RUN_OR_PICKLE
############################################################################
# run module and dump in pickle or retreive pickle without running module
def run_or_pickle(function_return_obj, function_name, arg1 = None, 
                  arg2 = None, arg3 = None, arg4 = None, arg5 = None):

  # make string with pickle name
  pickle_name = (function_return_obj + ".p")
  #logger.debug(pickle_name)

  # check if pickle file exists
  if os.path.isfile(pickle_name) == True:
    # retrieve pickle
    function_return_obj = pickle.load(
                          open((str(function_return_obj) +".p"),"rb"))
    logger.info('You have already analysed this bit, ' +
                'so we have pickled the ' + str(pickle_name))
  # otherwise we want to run the function
  else:
    # case 1: no arguments
    if (arg1 == None and arg2 == None and arg3 == None and arg4 == None and 
        arg5 == None):
      
      function_return_obj = function_name()

    # case 2: one argument
    elif arg2 == None:
      function_return_obj = function_name(arg1)

    # case 3: two arguments
    elif arg3 == None and arg4 == None and arg5 == None:
      function_return_obj = function_name(arg1, arg2)

    #case 4: all arguments
    else:
      function_return_obj = function_name(arg1, arg2, arg3, arg4, arg5)

    # dump result in pickle
    pickle.dump(function_return_obj, open(pickle_name, "wb"))

  # return what the function returned or the pickle
  return function_return_obj

############################################################################




############################################################################
### MAIN FUNCTION
############################################################################

def main():
  # greeting
  logger.info("Hi there, you are running drug_repo for drug repositioning!" +
              " Let's do some mapping.")


  # generate chembl dictionary
  chembl_dic = run_or_pickle("chembl_dic", process_chembl)

  # get list of uniprot ids from chembl_dic
  chembl_uniprot_list = run_or_pickle("chembl_uniprot_list", 
                                      flatten_values, chembl_dic)

  # generate drugbank_dictionary
  drugbank_dic = run_or_pickle("drugbank_dic", process_drugbank)

  # get list of uniprot ids from drugbank
  drugbank_uniprot_list = list(drugbank_dic)

  # merge lists and rm duplicates
  uniprot_list = run_or_pickle("uniprot_list", merge_lists, 
                              chembl_uniprot_list, drugbank_uniprot_list)

  #logger.debug(len(uniprot_list))
  
  ### OVERWRITE UNIPROT_LIST WITH MADE-UP LIST
  # overwrite the list with a small set ['B6DTB2', 'Q4JEY0','P11511']
  #['Q4JEY0', 'P68363', 'P10613', 'P18825', 'Q9UM73', 'E1FVX6']
  #uniprot_list = ['Q4JEY0', 'P68363', 'P10613','P18825']
  ###

  # run or pickle uniprot_to_arch to retrieve cath domain architectures
  cath_dic = run_or_pickle("cath_dic", uniprot_to_arch, uniprot_list, "cath")

  # generate list, flatten it and rm duplicates
  cath_list = run_or_pickle("cath_list", flatten_values, cath_dic)

  #logger.debug(len(cath_list))

  # run or pickle uniprot_to_arch to retrieve pfam domain architectures
  pfam_dic = run_or_pickle("pfam_dic", uniprot_to_arch, uniprot_list, "pfam")

  # generate list, flatten it and rm duplicates
  pfam_list = run_or_pickle("pfam_list", flatten_values, pfam_dic)


  # call archindex on cath values to find the ones from schisto
  uniprot_schisto_cath_dic = run_or_pickle(
      "uniprot_schisto_cath_dic", arch_to_uniprot, cath_list, "cath")
  
  # generate list, flatten it and rm duplicates
  uniprot_schisto_cath_list = run_or_pickle("uniprot_schisto_cath_list", 
                              flatten_values, uniprot_schisto_cath_dic)
  
  #logger.debug(len(uniprot_schisto_cath_list))

  # call archindex on pfam values to find ones from schisto
  uniprot_schisto_pfam_dic = run_or_pickle(
    "uniprot_schisto_pfam_dic", arch_to_uniprot, pfam_list, "pfam")
  
  # generate list, flatten it and rm duplicates
  uniprot_schisto_pfam_list = run_or_pickle("uniprot_schisto_pfam_list", 
                              flatten_values, uniprot_schisto_pfam_dic)
  #logger.debug(len(uniprot_schisto_pfam_list))


  # merge and rm duplicates
  # this is total list of unique schisto uniprot ids
  uniprot_schisto_list = run_or_pickle("uniprot_schisto_list", merge_lists, 
                        uniprot_schisto_cath_list, uniprot_schisto_pfam_list)

  logger.debug(len(uniprot_schisto_list))

  ### OVERWRITE UNIPROT_SCHISTO_LIST WITH MADE-UP LIST
  # this one is the 10 reviewd results ['P13566','Q9U8F1','P33676',
  #'P30114','Q26499','P16641','C4PZQ3','P37227','C4QCD2','Q5D8V5']
  #example of unreviewed entries: 'C7TY75', 'G4LXF4', 'C1L491'
  #uniprot_schisto_list = ['P13566','Q9U8F1','P33676',
  #'P30114','Q26499','P16641','C4PZQ3','P37227','C4QCD2','Q5D8V5']
  ###

  # filter list for only reviewed ones
  uniprot_schisto_filt = run_or_pickle("uniprot_schisto_filt", 
                                        expasy_filter,uniprot_schisto_list)

  logger.debug(uniprot_schisto_filt)

  # make list of cath ids that point to the uniprot_schisto_list

  #cath_final_dic = uniprot_to_arch(uniprot_schisto_filt, "cath")
  #logger.debug(cath_final_dic)
  #logger.debug(flatten_values(cath_final_dic))


  ### OVERWRITE FILTERED LIST
  #uniprot_schisto_filt = ['P33676']
  ###

  # generate big map for chembl drugs
  chembl_repo_dic = run_or_pickle(
    "chembl_repo_map", chembl_repo_map, chembl_dic, cath_dic, 
    uniprot_schisto_cath_dic, pfam_dic, uniprot_schisto_pfam_dic)

  logger.debug(chembl_repo_dic)

  # # loop over list of uniprot ids
  # for uniprot in uniprot_schisto_filt:
  #   # empty list for cath ids for each uniprot
  #   retro_cath = []
  #   # loop over cath dictionary
  #   for key in uniprot_schisto_cath_dic:
  #     # check if uniprot is in the values
  #     if uniprot in uniprot_schisto_cath_dic[key]:
  #       # append cath ids to the list
  #       retro_cath.append(key)
  #       retromap_dic[uniprot] = retro_cath



############################################################################




############################################################################
### call main function, prevent execution on import
############################################################################
if __name__ == "__main__":
  main()
############################################################################