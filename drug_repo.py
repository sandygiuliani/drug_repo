# Copyright 2014 Sandra Giuliani
# drug_repo.py

# A bioinformatics pipeline for drug repositioning for schistosomiasis

# Please dowload the full repository from github.com/sandygiuliani/drug_repo
# See README.md for more info




############################################################################
### IMPORT CONFIG
############################################################################
# import configuration file config.py, raise warning if absent

try:
  import config as c
# handle importerror if config.py is missing
except ImportError:
  logger.info('You are running drug_repo.py for drug repositioning.')
  logger.error('The configuration file config.py is missing' +
                    ' in the current directory!')
  logger.warning('The program is aborted.')
  # exit script
  sys.exit()
############################################################################




############################################################################
### LOGGER SET UP
############################################################################
# set up logger for logging messages and write to log log_drug_repo.log

import logging
# set up log file to write to, it will be overwritten every time ('w' mode)
# leave this level setting to DEBUG
logging.basicConfig(filename=c.log_name, filemode='w',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)
# leave this level setting to DEBUG
logger.setLevel(logging.INFO)
# create console handler
ch = logging.StreamHandler()
# CHANGE THIS TO TUNE LOGGING LEVEL from DEBUG/INFO/WARNING
ch.setLevel(logging.INFO)
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
### IMPORT PYTHON MODULES
############################################################################
# import python modules

# import os (old system) - deprecated, use subprocess instead
import os

# import os.path for checking if files exist
import os.path

# import pickle
#import pickle

# import regex
import re

# cpickle
import cPickle as pickle

# import csv for comma-sep files
import csv

# import subprocess for executing command line
import subprocess

# import itertools for flatten out lists
import itertools

# izip_longest to split dictionaries into chunks (part 8)
from itertools import izip_longest

# import other modules
import sys, re, string, fnmatch, shutil

# for http
from urllib2 import urlopen, HTTPError
#import urllib2

from urllib import urlretrieve
# linecache for jumping to specific lines
import linecache

# datetime
from datetime import datetime

# autovivification for creating nested dictionaries automatically
class AutoVivification(dict):
  """Implementation of perl's autovivification feature."""
  def __getitem__(self, item):
      try:
          return dict.__getitem__(self, item)
      except KeyError:
          value = self[item] = type(self)()
          return value

import gzip
import xml.dom.minidom
############################################################################




############################################################################
### IMPORT BIOPYTHON
############################################################################
# import biopython modules

# import Biopython Entrez
from Bio import Entrez
# tell NCBI who I am
Entrez.email = c.your_email

# list available databases
#handle = Entrez.einfo()
#logger.info(handle.read())

# import SeqIO
from Bio import SeqIO

# import expasy for accessing protein sequences
from Bio import ExPASy

# import swissprot for parsing swissprot plain text files
from Bio import SwissProt
############################################################################
  



############################################################################
### IMPORT MODELLER
############################################################################
# import modeller for homology modelling

if c.steps >= c.modeller_step:

  try:
    #load standard modeller classes
    from modeller import *
    # from modeller import automodel
    from modeller.automodel import *

  # handle error
  except ImportError:
    logger.info('You are running drug_repo.py for drug repositioning.')
    logger.error('MODELLER cannot be found!')
    logger.info('Please check it is properly installed' +
                ', or set a number of steps less than ' + 
                  str(c.modeller_step) +
                ' in the configuration file.')
    
    logger.warning('The program is aborted.')
    # exit script
    sys.exit()
############################################################################




############################################################################
### HEADER_COUNT
############################################################################
# find specific header in line of comma or tab (or other) separated
# file and return column number of the header

def header_count(line, separator, header):
  #set counter to 0
  col_number = 0
  # loop until word matches the header
  while not line.split(separator)[col_number].rstrip('\r\n') == header:
    col_number = col_number + 1
  # return column number
  return col_number
############################################################################



############################################################################
### FILE_TO_LINES
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
    logger.warning('The program is aborted.')
    # exit python
    sys.exit()
############################################################################




############################################################################
### SWAP_DIC HELPER
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
### PROCESS_CHEMBL
############################################################################
# process chembl input file, return dic of chembl ids vs uniprot ids

def process_chembl(input_file):
  # open chembldrugs.txt for reading
  lines = file_to_lines(input_file)
  # logger.info('The ChEMBL input file contains a total of '
  #             + str(len(lines)-1) + ' drugs.')

  ### ANALYSE HEADERS AND OBTAIN COLUMN NUMBERS
  logger.info(len(lines))
  # get the headers
  headers = lines[0]

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

  # logger.info('The number of drugs in phase 1 is: ' + str(phase1) +
  #             '; in phase 2 is: ' + str(phase2) + '; in phase 3 is: ' +
  #             str(phase3)+'; in phase 4 is: ' + str(phase4) +
  #             '; in unknown phase is: ' + str(phase_unknown) + '.')
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
    if (rowsplit2[col_phase] in c.chembl_phases):
      # append the stripped lines to the list
      stripped.append(lines[y])

  # empty list in which to store filtered chembl drug ids
  chembl_filt_list = []
  # loop over, note here there is no header
  for line in stripped:
    # tab separate
    rowsplit3 = line.split("\t")
    # check if they are of the desired drug type
    if rowsplit3[col_type] in c.chembl_mol_type:
      #logger.debug(rowsplit3[col_chemblid])

      # make list of chembl ids we are interested in
      chembl_filt_list.append(rowsplit3[col_chemblid])

  logger.info('We have filtered the entries in clinical phases ' +
              str(c.chembl_phases) + ' and of molecule type ' +
              str(c.chembl_mol_type) + ', to obtain ' +
              str(len(chembl_filt_list)) + ' drugs.')
  ###


  # CREATE DICTIONARY CHEMBL_TARGET_DRUG_DIC IN WHICH TO STORE
  # {CHEMBL TARGET IDS: (LIST OF CHEMBL DRUG IDS)}

  # open the drug targets chembl file and get lines
  drug_targ = file_to_lines(c.chembl_targets)
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
  # more than one chembl id can point to the same uniprot id!
  chembl_uniprot_map_dic = swap_dic(c.chembl_uniprot)
  #logger.debug(chembl_uniprot_map_dic)


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

  # return dictionary {uniprot1:(list of chembl ids)}
  return chembl_dic
############################################################################




############################################################################
### PROCESS_DRUGBANK
############################################################################
# process drugbank file, return a dictionary of drugs vs uniprot ids. Deal
# with duplicate values and append to list in the dictionary.

def process_drugbank(input_file):

  # open and read drug_bank input and count number
  lines = file_to_lines(input_file)
  # number of entries (remember there is header)
  # logger.info('The DrugBank input file contains ' + str(len(lines)-1) +
  #             ' entries.')

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
  for line in incsv:
    #list_check.append(line[col_uniprot])
    drug_string = line[col_drugbankids]
    # logger.info(drug_string)
    # list of drubbank ids
    drug_split = drug_string.split(';')
    # remove white spaces
    drug_split = [x.strip(' ') for x in drug_split]
    #ogger.debug(drug_split)
      # check if uniprot id is already in the dictionary
    if line[col_uniprot] in drugbank_dic:
      # append drug id value to the list in the dictionary
      # the '.extend' prevents the formation of a list of lists!
      drugbank_dic[line[col_uniprot]].extend(drug_split)
    else:
      # populate dictionary with the new entry
      drugbank_dic[line[col_uniprot]] = drug_split

  # confirm we are dealing with duplicates
  #logger.debug(len(list_check))
  #logger.debug(len(list(set(list_check))))


  # fast and silly way to generate the swap dictionary!

  drugbank_list = list(itertools.chain(*list(drugbank_dic.values())))
  drugbank_unique = list(set(drugbank_list))

  #logger.debug(len(drugbank_unique))



  drugbank_swap = {}

  # loop over list of unique drugbank ids
  for drug in drugbank_unique:
    #empty list for the uniprot values
    uniprot_list = []
    for uniprot in drugbank_dic:
      for db in drugbank_dic[uniprot]:
        if db == drug:
          # append uniprot value to the list
          uniprot_list.append(uniprot)

    # sanity check the list is not empty
    if uniprot_list:
    # the list is now ready, populate the dictionary
      drugbank_swap[drug] = uniprot_list



  #logger.debug(drugbank_swap)

  # return dictionary {uniprot1:(list of drugbank ids)}
  return drugbank_swap
############################################################################




############################################################################
### UNIPROT_TO_ARCH FUNCTION
############################################################################
# run archindex, return dictionary of CATH/pfam domain architecture vs
   # uniprot values

def uniprot_to_arch(uniprot_list,architecture):

  if architecture == "cath":
    flag = "-cath"

  elif architecture == "pfam":
    flag = ""

  # dictionary of uniprot ids and list of correposponding architectures
  arch_dic = {}
  numb_unipr = 0
  # loop over list of uniprot values
  for uniprot_id in uniprot_list:
    logger.info(uniprot_id)
    numb_unipr = numb_unipr + 1
    logger.info('Processing ' + str(numb_unipr))
    #list in which to store list of CATH domains for each entry
    architect_list = []
    # call archschema on the list
    subprocess.call(c.archindex_path + " -u " + str(uniprot_id) +
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
          # cath format eg '4.10.400.10'
          cath_format = re.compile('.*\..*\..*\..*')
          
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
              if cath_format.match(item):
                architect_list.append(item)

          # this is the case of just one entry, no undescores
          else:
            # check the format is CATH one
            if cath_format.match(line_nops):
              architect_list.append(line_nops)


        ### for pfam ###
        elif architecture == "pfam":
          # check if there are dots
          if "." in line_split[2]:
            dot_split = line_split[2].split(".")
            #logger.debug(undersc_split)

            for item in dot_split:

              architect_list.append(item)

          # this is the case of just one entry
          else:
            architect_list.append(line_split[2])

    logger.info(architect_list)
    # check list is not empty
    if architect_list:

      # eliminate duplicates within the list (this is for each entry!!)
      architect_list = list(set(architect_list))

      # populate the dictionary
      arch_dic[uniprot_id] = architect_list

  #logger.debug(cath_dic)
      # logger.info(arch_dic)
  # rm temp.txt in the end
  # this is the last temp file that overwrote the others
  subprocess.call("rm temp.txt", shell=True)

  # return dic
  return arch_dic
############################################################################




############################################################################
### ARCH_TO_UNIPROT FUNCTION
############################################################################
# run archindex, filter for TAXA, find uniprot ids, return arch vs uniprot dic

def arch_to_uniprot(arch_list,architecture):

  spam = 0

  if architecture == "cath":
    flag = "-cath"
  elif architecture == "pfam":
    flag = ""
  
  # regex starts with colon
  starts_colon = re.compile(':.*')


  # empty dictionary
  uniprot_dic = {}
  # loop over list of arch values
  for arch_id in arch_list:
    spam = spam + 1
    # logger.info(spam)
    # empty list in which to store uniprot values
    uniprot_list = []
    # iterate over the taxa code list (schisto species)
    for taxa_code in c.taxa:
      # call archschema on the list
      subprocess.call(c.archindex_path + " -p " + str(arch_id) +
                    " -maxa 1 -maxs 100 " + str(flag) + " -s " + taxa_code +
                    " > temp.txt", shell=True)
      # store lines
      lines = file_to_lines('temp.txt')
      #logger.info(lines)

      # get the uniprot values and append them to uniprot_list
      for i in range(len(lines)):
        # find line that starts with parent
        if lines[i][0:2] == ':A':
          
          counter = 0
          # while the line is not finishing line
          while not starts_colon.match(lines[i+2+counter].split("\t")[0]):
            # take the line' and split it
            line_split = lines[i+2+counter].split("\t")
            uniprot_list.append(line_split[0])
            counter = counter +1


    # check list is not empty
    if uniprot_list:
      # remove duplicates from list
      uniprot_list = list(set(uniprot_list))

      # populate the dictionary
      uniprot_dic[arch_id] = uniprot_list

  # rm temp.txt in the end
  # this is the last temp file that overwrote the others
  subprocess.call("rm temp.txt", shell=True)

  #return the dictionary
  return uniprot_dic
############################################################################




############################################################################
### MERGE_LISTS
############################################################################
# merge lists, remove duplicates and return merged list
def merge_lists(list1, list2):

  merged_list = list(set(list1)|set(list2))


  return merged_list
############################################################################




############################################################################
### EXPASY_FILTER
############################################################################
# take list of uniprot, call expasy and apply filter to find the reviewed
# entries (returns list)

  # entrez fetch does not seem to handle non-swissprot ids (TrEMBL)
  #prothandle = Entrez.efetch(
  #                  db="protein", id= 'P13566', rettype = "gb")
  #seq_record = SeqIO.read(prothandle, "gb")
  #logger.debug(seq_record)

def expasy_filter(uniprot_list, filter_type):
# filter_type == 'reviewed' if you want to keep only reviewed uniprots
# filter_type == 'pdb' if you want to keep only the ones with pdb structures

# alternative to swissport read, should also work
#record = SeqIO.read(handle, "swiss")


  # empty list to store filtered entries
  filtered_list = []

  # counter for obsolete/http err
  obsolete = 0

  for entry in uniprot_list:
    try:
        # get handle
        handle = ExPASy.get_sprot_raw(entry)
        # swissprot read
        record = SwissProt.read(handle)
        

        # FILTER according to presence pdb structure
        if filter_type == 'pdb':
          #logger.debug(record.cross_references)
          #so far the pdb does not exist
          pdb_exists = False 

          # get tuples with references (pdb, etc...)
          cross_ref = record.cross_references
          # loop over every tuple in the list
          for tup in record.cross_references:
            # if the first is pdb means there one or more pdb
            if tup[0] == "PDB":
              pdb_exists = True 

          if pdb_exists == True:
            # add entry to the list
            filtered_list.append(entry)


        # FILTER according to reviewed uniprot
        elif filter_type == 'reviewed':
          # check it is reviewed
          if record.data_class == 'Reviewed':
            # add reviewed entry to the list
            filtered_list.append(entry)


    except HTTPError:
      obsolete = obsolete + 1
      #logger.debug('uh-ho')


  logger.debug(obsolete)
  return filtered_list
############################################################################




############################################################################
### EXPASY_DIC
############################################################################
# take list of uniprot, run expasy, retrieve info and write dictionary
# eg {taxa1:[list of uniprot], taxa2:[list of uniprot]}

def expasy_dic(uniprot_list, filter_type):        

  # dictionary
  final_dic = {}

  for entry in uniprot_list:
    try:
        # get handle
        handle = ExPASy.get_sprot_raw(entry)
        # swissprot read
        record = SwissProt.read(handle)

        # FILTER according to reviewed uniprot
        if filter_type == 'taxa':
          # list of uniprots matching each
          taxa_id = record.taxonomy_id[0]
          # if it is already there
          if taxa_id in final_dic:
            # add entry to list
            final_dic[taxa_id].append(entry)

          else:
            entry_list = []
            entry_list.append(entry)
            # add entry (as a list) to dic
            final_dic[taxa_id] = entry_list
          
    except HTTPError:
      pass
  
  return final_dic

############################################################################




############################################################################
### FLATTEN_DIC
############################################################################
# list from dictionary's key OR values and flatten it 
# for the values, flatten list (from list of lists to simple list) 
# also eliminate duplicates and sort

def flatten_dic(dic, keys_or_val):
  # generate list from values
  if keys_or_val == "values":
    # generate list and flatten it
    flat_list = list(itertools.chain(*list(dic.values())))
  # generate list from keys
  elif keys_or_val == "keys":
    flat_list = dic.keys()

  # fourth level complex dic
  elif keys_or_val == "values_4":
    flat_list = []
    for a in dic:
      for b in dic[a]:
        for c in dic[a][b]:
          for d in dic[a][b][c]:
            flat_list.append(d)

  # rm duplicates
  flat_list = list(set(flat_list))

  #rm white spaces
  flat_list = [x.strip(' ') for x in flat_list]

  # sort
  flat_list.sort()

  return flat_list
############################################################################




############################################################################
### CHEMBL_REPO
############################################################################
# make big map for chembl ids, for both cath and pfam arch
# format {drug1:{drug_target1:{arch1:[list os schisto uniprot], arch2:[..]},
#                 drug_target2:{..}, drug2:{........}}}

# takes chembl to uniprot dic, uniprot to cath dic, cath to schisto dic,
# uniprot to pfam dic, pfam to schisto dic

def chembl_repo(chembl_dic, cath_dic, schisto_cath_dic,
                    pfam_dic, schisto_pfam_dic):

  # autovivification nested dictionary
  chembl_repo_map = AutoVivification()

  # loop over each drug in the dictionary
  for drug in chembl_dic:

    # dictionary for each drug
    each_drug_dic = {}

    # sanity check the uniprot list is not empty (it should not be empty)
    if chembl_dic[drug]:
      # loop on list of uniprot
      for target in chembl_dic[drug]:

        # make list in which to store the cath/pfam
        cath_pfam = []

        # append cath value, only if it exists in the dictionary
        if target in cath_dic:
          cath_pfam.extend(cath_dic[target])

        # append pfam value, only if it exists in the dictionary
        if target in pfam_dic:
          cath_pfam.extend(pfam_dic[target])

        #logger.debug(cath_pfam)
        # check the list is not empy
        if cath_pfam:
        # make dictionary
          each_drug_dic[target] = cath_pfam

      # finished loops for all targets, we have the dictionary
      # with {uniprot1:[list of arch],uniprot2:[list]}
      #logger.debug(each_drug_dic)

      #loop over each uniprot in new dictionary
      for targ in each_drug_dic:

        # dictionary for each target
        each_target_dic = {}
        # sanity check the arch list is not empty:
        if each_drug_dic[targ]:
          # loop over list of architectures
          for arch in each_drug_dic[targ]:
            # empty list for schisto
            schisto_list = []
            # check it exists in the cath dictionary
            if arch in schisto_cath_dic:
              schisto_list.extend(schisto_cath_dic[arch])

            if arch in schisto_pfam_dic:
              schisto_list.extend(schisto_pfam_dic[arch])



            # check list is not empty
            if schisto_list:
              #logger.debug(schisto_list)
              # add list as values in dictionary
              each_target_dic[arch] = schisto_list

              #logger.debug(each_target_dic)


              chembl_repo_map[drug][targ][arch] = each_target_dic[arch]


  return chembl_repo_map
############################################################################




############################################################################
### DRUGBANK_REPO
############################################################################
# make big map for drugbank ids, for both cath and pfam arch
# format {drug1:{drug_target1:{arch1:[list os schisto uniprot], arch2:[..]},
#                 drug_target2:{..}, drug2:{........}}}

# takes chembl to uniprot dic, uniprot to cath dic, cath to schisto dic,
# uniprot to pfam dic, pfam to schisto dic

def drugbank_repo(drugbank_dic, cath_dic, schisto_cath_dic,
                    pfam_dic, schisto_pfam_dic):

  # autovivification nested dictionary
  drugbank_repo_map = AutoVivification()

  #below copied from chembl_repo_map

  # loop over each drug in the dictionary
  for drug in drugbank_dic:

    # dictionary for each drug
    each_drug_dic = {}

    # sanity check the uniprot list is not empty (it should not be empty)
    if drugbank_dic[drug]:
      # loop on list of uniprot
      for target in drugbank_dic[drug]:

        # make list in which to store the cath/pfam
        cath_pfam = []

        # append cath value, only if it exists in the dictionary
        if target in cath_dic:
          cath_pfam.extend(cath_dic[target])

        # append pfam value, only if it exists in the dictionary
        if target in pfam_dic:
          cath_pfam.extend(pfam_dic[target])

        #logger.debug(cath_pfam)
        # check the list is not empy
        if cath_pfam:
        # make dictionary
          each_drug_dic[target] = cath_pfam

      # finished loops for all targets, we have the dictionary
      # with {uniprot1:[list of arch],uniprot2:[list]}
      #logger.debug(each_drug_dic)

      #loop over each uniprot in new dictionary
      for targ in each_drug_dic:

        # dictionary for each target
        each_target_dic = {}
        # sanity check the arch list is not empty:
        if each_drug_dic[targ]:
          # loop over list of architectures
          for arch in each_drug_dic[targ]:
            # empty list for schisto
            schisto_list = []
            # check it exists in the cath dictionary
            if arch in schisto_cath_dic:
              schisto_list.extend(schisto_cath_dic[arch])

            if arch in schisto_pfam_dic:
              schisto_list.extend(schisto_pfam_dic[arch])

            # check list is not empty
            if schisto_list:
              #logger.debug(schisto_list)
              # add list as values in dictionary
              each_target_dic[arch] = schisto_list

              #logger.debug(each_target_dic)

              drugbank_repo_map[drug][targ][arch] = each_target_dic[arch]

  #logger.debug(drugbank_repo_map)
  return drugbank_repo_map
############################################################################




############################################################################
### FILT_SCHISTO_MAP
############################################################################
# filter the chembl_repo_map for the filtered schisto uniprots we are
# interested in

def filt_schisto_map(chembl_repo_map, schisto_filt):
  # empty dictionary

  # autovivification nested dictionary
  schisto_filt_map = AutoVivification()

  for drug in chembl_repo_map:
    #schisto_filt_map[drug] = {}
    for target in chembl_repo_map[drug]:
      #schisto_filt_map[drug][target] = {}
      for arch in chembl_repo_map[drug][target]:
        #logger.debug(arch

        schisto_list = []

        for schisto in chembl_repo_map[drug][target][arch]:
          if schisto in schisto_filt:
            schisto_list.append(schisto)
            if schisto_list:
              schisto_filt_map[drug][target][arch] = schisto_list

  return schisto_filt_map
############################################################################




############################################################################
### CSV_TO_DIC
############################################################################
# take csv file such as uniprot_pdb, return dictionary with first column
# values as keys and ';'-separated values as list

def csv_to_dic(csv_file):


  # empty dictionary
  csv_dic = {}

  # open and read drug_bank input and count number
  lines = file_to_lines(csv_file)

  # read csv lines with the csv reader - deals with quotation marks etc..
  incsv = csv.reader(lines)
  # skip the first line with the headers
  next(incsv)

  # loop over each line
  for line in incsv:
    # empty list to store the pdb ids
    pdb_list = []
    pdb_strip = line[1].split(';')
    for pdb in pdb_strip:
      pdb_list.append(pdb)
  
   #   logger.debug(line)
    csv_dic[line[0]] = pdb_list

  #logger.debug(len(csv_dic))

  return csv_dic

############################################################################




############################################################################
### CSV_TO_LST
############################################################################
# take csv file and return single list of the last line
# it does not properly skip the comment lines! 
# it works but better ADD WRAPPER TO PROPERLY SKIP THE COMMENT LINE!

def csv_to_lst(csv_file):
  # regular expression for string containing at least one '#'
  #contains_comment = re.compile('#.*')

  # # open and read file
  # lines = file_to_lines(csv_file)
  # # count comment lines
  # comment_count = 0
  # # loop over each line
  # for line in lines:
  #   logger.debug(line)
  #   # skip comment lines
  #   if CONTAINS_COMMENT.match(line):
  #     comment_count = comment_count + 1
  # logger.debug(comment_count)


  with open(csv_file) as f:

    #,quoting=csv.QUOTE_MINIMAL
    # get rows with csv reader
    # the quotechar is for ignoring the single quotes I have in the list
    # otherwise you end up with double quotes in the list
    rows = csv.reader(f, delimiter=',', quotechar='\'', skipinitialspace=True)
    for row in rows:
      # overwrites until it reaches last row
      csv_list = row
      #logger.debug(row)
    
    # order list and rm duplicates
    csv_list = list(set(csv_list))
    csv_list.sort()

  # logger.debug(csv_list)
  # return the list
  return csv_list

############################################################################




############################################################################
### FILTER_DIC_FROM_LIST
############################################################################
# take dictionary and return filtered dictionary with only entries that are
# in the given list

def filter_dic_from_list(dictionary,filt_list):
  filtered_dic = {}

  for item in dictionary:
    if item in filt_list:
      filtered_dic[item] = dictionary[item]

  return filtered_dic

############################################################################




############################################################################
### LIST_SECOND_LEVEL_DIC
############################################################################
# take dictionary of dictionaries and return unique list of the second level
# also for values of a simple dic

def list_second_level_dic(dictionary):
  second_list = []

  for item in dictionary:
    for sec in dictionary[item]:
      # strip white spaces
      sec.strip()
      second_list.append(sec)

  second_list = list(set(second_list))

  return second_list

############################################################################




############################################################################
### LST_DIC
############################################################################
# read lst file and creat dictionary with first column as key and rest as list
def lst_dic(lst_file):
  lines = file_to_lines(lst_file)
  lst_dictionary = {}
  # iterate over lines
  for i in range(0,len(lines)):
    #logger.debug(line)
    # split tab
    splitline = lines[i].split(":")
    #logger.debug(splitline)
    
    # empty list for the het groups
    het_list = []

    # remove (right strip) carriage return and the final ';'
    het_gr = splitline[1].rstrip('\r\n').rstrip(';')
    #het_gr = het_gr.rstrip(' ')
    #logger.debug(het_gr)
    het_split = het_gr.split(";")
    #logger.debug(het_split)
    for het_group in het_split:
      het_list.append(het_group.strip())
    #logger.debug(het_list)
    # populate dictionary, stripping the carriage return
    lst_dictionary[splitline[0].strip()] = het_list
  
  #logger.debug(lst_dictionary)
  return lst_dictionary
############################################################################




############################################################################
### EXCLUDE_VALUES_FROM_DIC
############################################################################
# take dictionary and return filtered dictionary excluding values that do
# come up in list
# flag exclude - keep values that are not in the list
# flag include - only keep values that are in the list
# flag nomatch - keep values that do not match the expression

def exclude_values_from_dic(dictionary, filt_list, flag):
  filtered_dic = {}

  for item in dictionary:
    pass_list = []
    for value in dictionary[item]:
      # if we want to exclude the ones that are in the list
      if flag == 'exclude':
        if value not in filt_list:
          pass_list.append(value)
      # if we want only the ones that are in the list
      elif flag == 'include':
        if value in filt_list:
          pass_list.append(value)

      elif flag == 'nomatch':
        if not filt_list.match(value):
          pass_list.append(value)
    # check the list is not empty
    if pass_list:
      # popolate filtered dic
      filtered_dic[item] = pass_list

  return filtered_dic

############################################################################




############################################################################
### TXT_TO_DIC FUNCTION
############################################################################
# take txt (comma-sep, with headers) and create dictionary of header 1 
# vs header 2

def txt_to_dic(input_file, header1, header2):
  # open chembldrugs.txt for reading
  lines = file_to_lines(input_file)
  # get the headers
  headers = lines[0]
  # logger.debug(headers)

  # empty dic
  dic = {}

  # number of header 1 nad 2
  h1 = header_count(headers, "\t", header1)
  h2 = header_count(headers, "\t", header2)

  # populate dictionary
  for i in range(1,len(lines)):
    dic[lines[i].split("\t")[h1]] =  lines[i].split("\t")[h2].rstrip('\r\n')
  

  return dic

############################################################################




############################################################################
### SMI_TO_DIC FUNCTION
############################################################################
# take smi (comma-sep, without headers) and create dictionary of header 1 
# (column number n1) vs header 2 (column number n2)

def smi_to_dic(input_file, n1, n2):
  # open for reading
  lines = file_to_lines(input_file)
  #logger.debug(lines[18482])

  # empty dic
  dic = {}

  # populate dictionary
  for i in range(0,len(lines)):
    # check line is not empty string
    split = lines[i].rstrip('\r\n').split("\t")
    
    # check the number of items is the same as the number of items in the
    # first row of the file
    if len(split) == len(lines[0].split("\t")): 
    #or len(split) == (len(lines[0].split("\t"))-1):

      if split[n1].rstrip('\r\n') != '' and split[n2].rstrip('\r\n') != '':
        #logger.debug(split)
        dic[split[n1]] =  split[n2]

  # lenght total file (nb there are new lines!)      
  # logger.debug(len(lines))
  # # lenght dictionary created
  # logger.debug(len(dic))

  return dic

############################################################################




############################################################################
### DIC_TO_TXT
############################################################################
# take dic and write to file as tab-sep file

def dic_to_txt(dic, file_name):

  #try remove file if it exists already
  if os.path.isfile(file_name) == True:
    os.remove(file_name)


  with open(file_name, 'a') as myfile:

    for item in dic:
      line = str(dic[item])+ "\t" + str(item) + "\r\n"
      # logger.debug(line)
      myfile.write(line)
      #myfile.write('\n')
      # don't forget you will have new line at the end!

############################################################################




############################################################################
### BABEL_SMI_TO_SDF
############################################################################
# call open babel to convert smiles into sdf
# add --gen3d for coordinates + " --gen3d"
# -d to delete hydrogens

def babel_smi_to_sdf(input_file, output_file):
  # check if the output exists already
  if os.path.isfile(output_file) == True:
    logger.info('The converted file ' + str(output_file) + 
                ' exists already. If you wish to run the conversion again,' +
                ' delete it.')

  else:
   subprocess.call("babel" + " -ismi " + str(input_file) +
                  " -osdf " + str(output_file) + " -d --gen2d", shell=True)

  # nb this module returns None! it just writes to output.sdf

############################################################################




############################################################################
### SDF_TO_DIC
############################################################################
# takes sdf file and converts in dictionary with key versus value
# nb key has to be unique (eg drugbank identifier), 
# otherwise will be overwritten
def sdf_to_dic(sdf_file, key, value):
  with open(sdf_file) as f:
    line_count = 0
    # empty dic
    dic_from_sdf = {}

    for line in f:
      line_count = line_count + 1
      if key in line:
        db_id = linecache.getline(sdf_file, (line_count+1))
        # strip newline
        db_id = db_id.rstrip('\n')

      if value in line:
        smi_id = linecache.getline(sdf_file, (line_count+1))
        # strip newline
        smi_id = smi_id.rstrip('\n')
        dic_from_sdf[db_id] = smi_id

  #logger.info(len(dic_from_sdf))

  return dic_from_sdf

############################################################################




############################################################################
###  MERGE_DIC
############################################################################
def merge_dic(dic1,dic2,dic3):

  dic_merge = {}
  # each drug
  for item in dic1:
    #list for pdbs
    list2 = []
    # list for cc
    list3 = []

    # each uniprot
    for value in dic1[item]:
      if value in dic2:
        # append pdb values to list
        list2.extend(dic2[value])
        #logger.info(value)
        #logger.info(dic2[value])

    #logger.info(list2)
    # check list is not empty
    if list2:
      #rm duplicates
      list2 = list(set(list2))
      #for each pdb in this list
      for thing in list2:
        # check the pdb is in the dic3
        if thing in dic3:
          # append to list 3
          list3.extend(dic3[thing])


    # before moving to next drug, populate dic
    if list3:
      # rm duplicates
      list3 = list(set(list3))
      item = item.strip(' ')
      dic_merge[item] = list3
  #logger.info(len(dic_merge))

  return dic_merge
############################################################################




############################################################################
### RUN_SMSD
############################################################################
# call SMSD

def run_smsd(query, target, flag, threshold, dic_map=None):
  # flag = 'pair' for pairwise comparison, query and target are dictionaries
  # ids: smiles
  # flag = 'batch' for batch comparison, query is dictionary,
  # target is sdf file

  # run SMSD
  # can also include cwd=SMSD_PATH together with shell=True in the command
  # for instance:
  # sh ./SMSD.sh -Q SMI -q \"CCCCC\" -T SMI -t \"CCCN\" -O SMI -o test.smi
  # -r remove hydrogens
  # -m produce mapping output (molDescriptors.out, mcs.out and .mol files)
  # -b match bond type (bond sensitive, faster run!)
  # -z match rings (this is needed otherwise very slow)
  # -x match atom type... (?)
  # -g for png image
  # recommended options are -r -z -b

  # increase Java max heap size
  # subprocess.Popen(["export JVM_ARGS=\"-Xms1024m -Xmx1024m\""], shell=True)
  # subprocess.Popen(["java -Xmx1024m ..."], shell=True)

  # current directory
  initial_dir = os.getcwd()
  #logger.debug(initial_dir)

  if flag == 'pair_2dic':
    # navigate to SMSD dir
    os.chdir(c.smsd_path)

    output = open('smsd_run_pair_2dic.txt', 'w')

    # dictionary drugs: list of cc that match
    drug_cc_dic = {}

    drug1 = {}
    drug09 = {}
    drug08 = {}
    drug07 = {}

    #loop over drug in the map
    for drug in dic_map:
      #logger.info(drug)
      # list of cc that match each drug
      match_cc_list = []
      sim1 = []
      sim09 = []
      sim08 = []
      sim07 = []
     
      # check drug is in dic
      if drug in query:
        # get the drug smiles
        drug_smi = query[drug]

        
        # each cc we have in the list
        for cc in dic_map[drug]:
          #logger.info(cc)
          # check if cc is in the dic
          if cc in target:
            # get the smiles
            cc_smi = target[cc]

            # handles errors..
            try:
              subprocess.call("sh SMSD.sh -Q SMI -q \"" + str(drug_smi) + \
                            "\" -T SMI -t \"" + str(cc_smi) + 
                            "\" -r -m -z -b", shell=True)
              # get list from lines in molDescriptors output
              mol_des = file_to_lines("molDescriptors.out")
              #logger.info(mol_des)
        

              # make string out of list
              mol_str = ''.join(mol_des)
              # logger.info(mol_str)
              # sanity check string not empty
              if mol_str != '':

                # tanimoto similarity string or other string
                str_match = 'Tanimoto (Sim.)= '
                # find the string
                str_numb = mol_str.find(str_match)
                sim_index = (str_numb + len(str_match))

                sim_string = mol_str[sim_index:(sim_index+3)]
                # get similarity number and convert to float
                similarity = float(sim_string)

                if similarity >= float(threshold):
                  # append to list
                  match_cc_list.append(cc)

                if similarity >= 1.0:
                  sim1.append(cc)

                if similarity >= 0.9:
                  sim09.append(cc)

                if similarity >= 0.8:
                  sim08.append(cc)

                if similarity >= 0.7:
                  sim07.append(cc)
            except KeyboardInterrupt:
              logger.warning('You are terminating the script!')
              sys.exit()

            except:
              logger.warning('We have skipped one comparison!')
              pass


      # check the list is not empty
      if match_cc_list:
        # add list to dictionary
        drug_cc_dic[drug] = match_cc_list
        #logger.info(drug_cc_dic)

      if sim1:
        drug1[drug] = sim1
      if sim09:
        drug09[drug] = sim09
      if sim08:
        drug08[drug] = sim08
      if sim07:
        drug07[drug] = sim07  
      
    # write to output
    output.write("# similarity 1.0\n")
    output.write(str(drug1)+ "\n")
    output.write("# similarity > 0.9\n")
    output.write(str(drug09)+ "\n")
    output.write("# similarity > 0.8\n")
    output.write(str(drug08)+ "\n")
    output.write("# similarity > 0.7\n")
    output.write(str(drug07)+ "\n")
  
  # close the file
  output.close()


  # ###
  # # pairwise comparison on smile strings
  # ###
  # elif flag == 'pair':
  #   # navigate to SMSD dir
  #   os.chdir(SMSD_PATH)
  #   #logger.info('this is pairwise')
  #   # dictionary drugs: list of cc that match
  #   drug_cc_dic = {}

  #   for drug in query:
  #     drug_smi = query[drug]
  #     #logger.info(drug)
  #     # list of cc that match each drug - with similarity above threshold
  #     match_cc_list = []
  #     for cc in target:
  #       cc_smi = target[cc]
  #       # logger.info(cc)

  #       subprocess.call("sh SMSD.sh -Q SMI -q \"" + str(drug_smi) + \
  #                     "\" -T SMI -t \"" + str(cc_smi) + "\" -r -m -z -b", 
  #                     shell=True)

  #       # get list from lines in molDescriptors output
  #       mol_des = file_to_lines("molDescriptors.out")
  #       #logger.info(mol_des)
  #       # make string out of list
  #       mol_str = ''.join(mol_des)
  #       # logger.info(mol_str)
  #       # tanimoto similarity string or other string
  #       str_match = 'Tanimoto (Sim.)= '
  #       # find the string
  #       str_numb = mol_str.find(str_match)
  #       sim_index = (str_numb + len(str_match))
  #       # get similarity number and convert to float
  #       similarity = float(mol_str[sim_index:(sim_index+3)])

  #       if similarity >= float(threshold):
  #         match_cc_list.append(cc)

  #     # check the list is not empty
  #     if match_cc_list:
  #       # add list to dictionary
  #       drug_cc_dic[drug] = match_cc_list
    

  # ###
  # # batch processing on sdf
  # ###
  # elif flag == 'batch':
  #   logger.info('this is batch')
  #   # dictionary drugs: list of cc that match
  #   drug_cc_dic = {}
  #   # check if the target file is in the current directory
  #   if os.path.isfile(target) == False:
  #     logger.info('uh-oh')
  #     logger.error('The file ' + target + ' cannot be found' +
  #                  ' in the current directory!')
  #     # warning
  #     logger.warning('The program is aborted.')
  #     # exit python
  #     sys.exit()
  #     # logger.debug("cp " + target + " " + SMSD_PATH + "/" + target)

  #   else:
  #     logger.info('alright, the file is here')
  #     # copy the file to the SMSD directory
  #     subprocess.call("cp " + target + " " + SMSD_PATH + "/" + target, 
  #      shell=True)

  #     # move to SMSD directory
  #     os.chdir(SMSD_PATH)
  #     directory = os.getcwd()
  #     logger.debug(directory)

  #     for drug in query:
  #       # list of cc that match each drug - with similarity above threshold
  #       match_cc_list = []
  #       drug_smi = query[drug]
  #       logger.info(drug)

  #       subprocess.call("sh SMSD.sh -Q SMI -q \"" + str(drug_smi) + 
  #                     "\" -T SDF -t " + str(target) + " -m -r -z -b", 
  #                     shell=True)

  #       # get the ones above threshold

  #       # add to dic

  #       # check the list is not empty
  #       if match_cc_list:
  #         # add list to dictionary
  #         drug_cc_dic[drug] = match_cc_list


  # move back to working directory!
  os.chdir(initial_dir)
  directory = os.getcwd()
  #logger.debug(directory)
  
  # return the dictionary
  return drug_cc_dic

############################################################################




############################################################################
### MV_FILE
############################################################################
#take file, make copy and move it to current directory
# check if file is already there first!

def mv_file(origin, filename, new_name):

  # if the file is not there already
  if os.path.isfile(new_name) == False:
    current_dir = os.getcwd()

    os.chdir(origin)
    subprocess.call("cp " + filename + " " + current_dir + "/" + new_name, 
         shell=True)

    os.chdir(current_dir)


############################################################################




############################################################################
### FILTER_TXT
############################################################################
# write output file after filtering input_file according to filt_list

def filter_txt(input_file, output_file, header_name, filt_list):
  # open chembldrugs.txt for reading
  lines = file_to_lines(input_file)

  # get the headers
  headers = lines[0]
  # find header name
  col_header = header_count(headers, "\t", header_name)

  out = open(output_file, 'w')
  out.write(headers)
  #out.write('\n')
  for i in range(1,len(lines)):
    rowsplit = lines[i].split("\t")

    # get the chembl drug and target id values for the row
    col = rowsplit[col_header]

    # only proceed if the target id is not an empty field!
    if col  != "":

      # check if the molecule chembl id is one of the drugs we want
      if col in filt_list:
        out.write(lines[i])
        #out.write('\n')

  # close output file
  out.close()

############################################################################




############################################################################
### RUN_MODELLER
############################################################################
# run modeller, produce homology model
# this module would return None!

# alignment file name, code of template, code of sequence
def run_modeller(no_models, alnfile, knowns, sequence):

  # request minimal/verbose/none output
  log.none()    

  # create a new MODELLER environment to build this model in
  env = environ()  

  # directories for input atom files
  env.io.atom_files_directory = ['.']

  a = automodel(env, alnfile, knowns, sequence, 
                assess_methods=(assess.DOPE, assess.GA341))
  

  # (determines how many models to calculate)
  a.starting_model= 1                
  # index of last model
  a.ending_model  = no_models            

  # make the model
  a.make()

  # get list of successfully but models
  ok_models = filter(lambda x: x['failure'] is None, a.outputs)

  # rank models by DOPE score
  key = 'DOPE score'
  ok_models.sort(lambda a,b: cmp(a[key], b[key]))                        

  # return 'foo' 
  return ok_models  

############################################################################




############################################################################
### UNIPROT_TO_FASTA
############################################################################
# takes uniprot list and writes to file all fasta sequences
# returns foo
def uniprot_to_fasta(uniprot_list):

  # make string for file name
  file_name = str('dr_' + uniprot_list[0] + '_align.fasta')

  logger.info(file_name)
  
  with open(str(file_name), 'w') as f:
    # for uniprot in uniprot_list:
    for entry in uniprot_list:
      try:
          # get handle
          handle = ExPASy.get_sprot_raw(entry)
          # swissprot read
          record = SwissProt.read(handle)

          # list of uniprots matching each
          seq = record.sequence
          # add first line
          f.write('>' + entry + '\n')
          # add sequence
          f.write(seq + '\n')
          #logger.info(seq)
            
      except HTTPError:
        logger.error('We cannot access ExPASy!')
        logger.warning('We are exiting the program.')
        sys.exit()

  # simply return the file name
  return file_name


############################################################################


############################################################################
### RUN_TCOFFEE
############################################################################
# take fasta (fasta_file), run tcoffee to produce multisequence alignment

def run_tcoffee(fasta_file):
  # run t-coffeee
  # st out redirected to log file
  # - output=clustalw (default)
  # clustalw_aln, clustalw : ClustalW format.
  # gcg, msf_aln       : MSF alignment.
  # pir_aln  : pir alignment.
  # fasta_aln: fasta alignment.
  # phylip: Phylip format.
  # pir_seq: pir sequences (no gap).
  # fasta_seq: fasta sequences (no gap
  # quiet
  # -quiet=t_coffee.log
  try:
    subprocess.call("t_coffee " + str(fasta_file) + 
                    " -email=" + str(c.your_email) + 
                    " -quiet=" + str(c.t_coffee), shell=True)
    # get list from lines in molDescriptors output


  except KeyboardInterrupt:
    logger.warning('You are terminating the script!')
    sys.exit()

  except:
    logger.warning('T-coffee is failing, we are exiting the script!')
    sys.exit()

############################################################################




############################################################################
### RUN_OR_PICKLE
############################################################################
# run module and dump in pickle or retrieve pickle without running module

def run_or_pickle(function_return_obj, function_name, arg1 = None,
                  arg2 = None, arg3 = None, arg4 = None, arg5 = None):

  # make string with pickle name
  pickle_name = (function_return_obj + ".p")

  # check if pickle file exists
  if os.path.isfile(pickle_name) == True:
    # retrieve pickle
    function_return_obj = pickle.load(
                          open((str(function_return_obj) +".p"),"rb"))
    logger.debug('We have pickled ' + str(pickle_name))
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

    # case 4: three arguments
    elif arg4 == None and arg5 == None:
      function_return_obj = function_name(arg1, arg2, arg3)

    # case 5: four arguments
    elif arg5 == None:
      function_return_obj = function_name(arg1, arg2, arg3, arg4)

    #case 6: all arguments
    else:
      function_return_obj = function_name(arg1, arg2, arg3, arg4, arg5)

    # dump result in pickle
    pickle.dump(function_return_obj, open(pickle_name, "wb"))

  # return what the function returned or the pickle
  return function_return_obj

############################################################################




############################################################################
### TAXA_TO_SPECIES
############################################################################
# return list of species from taxa ids, to be used for logging purposes

def taxa_to_species(taxa_list, species_map):

  species_lst = []

  # string to be used in the loggers
  # species = ''

  # regex starts with 'N='
  starts_n = re.compile('N=.*')  

  for tax in taxa_list:

    bits_list = []
    with open(species_map) as f:

      for line in f:
        #split
        split_line = line.split(' ')
        if split_line[0].strip(' ') == tax:
          #print split_line
          for i in range(1,len(split_line)):
            #print split_line[i]
            if starts_n.match(split_line[i]):
              # print split_line[i]
              match = i
              for j in range(match,len(split_line)):
                bits_list.append(split_line[j].rstrip('\n').lstrip('N='))
              # print split_line[i]

    species_string = " ".join(bits_list)
    # print species_string
    species_lst.append(species_string)


  # return string of species
  return species_lst
############################################################################




############################################################################
### LST_TO_STRING
############################################################################
  
def lst_to_string(input_list):

  # string
  final_string = ''


  # string 'species' lists the species
  if len(input_list) == 1:
    final_string = str(input_list[0])

  if len(input_list) == 2:
    final_string = (input_list[0] + ' and ' + input_list[1])

  if len(input_list) > 2:
    for i in range(0,(len(input_list)-2)):
      final_string = (final_string + input_list[i] + ', ') 

    final_string = (final_string + input_list[len(input_list)-2] + ' and ')
    final_string = (final_string + input_list[len(input_list)-1])


  return final_string

############################################################################




############################################################################
### RES_NUMB_MAP
############################################################################
# return 2 residue numbers maps
# umap = uniprot:pdb /  pmap = pdb:uniprot
# module modified from 'siftPasser' (Nick Furnham)

def res_numb_map(doc):
  umap = {}
  pmap = {}
  c = {}
  d = {}
  pres = None
  ures = None
  #Data struc: {chain:{uni_res:pdb_res}}
  for e in doc.childNodes[0].childNodes:
    if e.nodeType == e.ELEMENT_NODE and e.localName == "entity":
      chain = ""
      #chain = str(e._attrs["entityId"].value)
      for f in e.childNodes:
        if f.nodeType == f.ELEMENT_NODE and f.localName == "segment":
          for g in f.childNodes:
            if (g.nodeType == g.ELEMENT_NODE and 
                g.localName == "listResidue"):
              for h in g.childNodes:
                if h.localName == "residue":
                  for i in h.childNodes:
                    if (i.localName == "crossRefDb" 
                        and i._attrs["dbSource"].value == "PDB"):
                      try:
                        pres = int(i._attrs["dbResNum"].value)
                        chain = str(i._attrs["dbChainId"].value)
                      except:
                        pres = int(i._attrs["dbResNum"].value[:-1])
                        chain = str(i._attrs["dbChainId"].value)
                    if (i.localName == "crossRefDb" and 
                        i._attrs["dbSource"].value == "UniProt"):
                      ures = int(i._attrs["dbResNum"].value)
                  if ures != None and pres != None:
                    c[ures]=pres
                    d[pres]=ures
                    ures = None
                    pres = None
          
      umap[chain] = c
      pmap[chain] = d
      c = {}
      d = {}

  return umap, pmap
############################################################################




############################################################################
### DRUG_TARG_RES_FILTER
############################################################################
# from drug mapping file get list or residue numbers (uniprot numbering)
# that interact with the drug
# drug_het_map drug:target:pdb:het

def drug_targ_res_filter(drug_het_map, drug_arch_target):
 
  # get pfam
  # these will be pdb numbering!!
  pfam_lines = file_to_lines(c.pdb_to_pfam)
  
  # get CATH
  # this will be uniprot numbering!!
  cath_lines = file_to_lines(c.uniprot_cath)

  # overwrite drug_het_map for testing!!
  # drug_het_map = 
  # {'CHEMBL1560':{'P12821': {'4c2p': ['X8Z'], '1uzf': ['MCO']}}}
  
  # drug:targt:cath:schistotarg
  drug_filt = AutoVivification()

  # logger.info(drug_het_map)
  for drug in drug_het_map:
    logger.info(drug)
    # list in which to store res
    res_list = []
    
    for target in drug_het_map[drug]:
      logger.info(target)

      # list of good domains! to save because there is
      # at least one of the pdbs that has that domain interacting with drug!
      # keep adding for each pdb
      good_dom = []

      for pdb in drug_het_map[drug][target]:
        logger.info('we are looking at pdb ' + str(pdb))
        
        pdb_upper = pdb.upper()
        
    
        # logger.info(pdb_upper)
        # dict for het to chain ID to residue
        het_ch_res = AutoVivification()
        
        ####################################
        # retreive SIFTS for uniprot-pdb res mapping
        ftp_url = ("ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/" + 
                      pdb + ".xml.gz")
        # gzip will not open ftp, urlretrieve it's needed
        f = gzip.open(urlretrieve(ftp_url)[0])
        doc = xml.dom.minidom.parse(f)
        # logger.info(doc)
        f.close()

        # get the two dics
        uni_pdb, pdb_uni = res_numb_map(doc)
        # logger.info(uni_pdb)
        ####################################


        #call pdbsum
        # logger.info(pdb)
        pdb_m = str(pdb[1]+pdb[2])
        # logger.info(pdb_m)
        # logger.info(drug_het_map[drug][target][pdb])
        psum = urlopen('http://www.ebi.ac.uk/thornton-srv/databases/PDBsum/' +
                       pdb_m + '/' + pdb + '/grow.out')
        
        psum_read = psum.readlines()
        # logger.info(psum_read)
        

        for i in range(1,len(psum_read)):
          # logger.info(psum_read[i])
          # split tab
          splitline = psum_read[i].split(" ")
          splitline = filter(None, splitline)
          het_name = splitline[11]
          chain_name = splitline[3]
          resnum = splitline[4] 
          
          # line with het group!
          # logger.info(drug_het_map[drug][target][pdb])
          # check the het gropu is one of the ones we want
          if het_name in drug_het_map[drug][target][pdb]:

            
            # check if it is already in there
            if het_name in het_ch_res:

              # check if chain name is already there
              if chain_name in het_ch_res[het_name]:
                # check if residue number is already there
                if resnum in het_ch_res[het_name][chain_name]:
                  # do nothing, it is already listed
                  pass
                # the chain is there, we need to add the res
                else:
                  het_ch_res[het_name][chain_name].append(resnum)

              # chain is not there
              else:
                res_list = []
                # res_list.append(resnum)
                het_ch_res[het_name][chain_name] = res_list
            
            # het is not in the dic, add res as list
            else:
              res_list = []
              # res_list.append(resnum)
              het_ch_res[het_name][chain_name] = res_list

          # het name is not of interested, skip the line
          else:
            pass


        pfam_numb = AutoVivification()
        #pfam - pdb numbering!!!
        for pfam_line in pfam_lines:

          line_split = pfam_line.split("\t")
          if line_split[0] == pdb_upper:
            # logger.info(line_split)
            # get the pfam id
            pfam_id = line_split[4].split(".")[0]
            p_chain = line_split[1]
            # logger.info(line_split[2])
            if line_split[2].isdigit():

              p_start = int(line_split[2])

              if line_split[3].isdigit():
                p_end = int(line_split[3])
                pfam_numb[pfam_id][p_chain] = (p_start, p_end)

              else:
                logger.info(str(line_split[3]) + '!!!!')
            
            else:
              logger.info(str(line_split[2]) + '!!!!')
            # p_end = int(line_split[3])
            # # logger.info(pfam_id)
            # pfam_numb[pfam_id][p_chain] = (p_start, p_end)

        # logger.info(pfam_numb)


        # get CATH - uniprot numbering!! -> is converted into pdb
        cath_numb = AutoVivification()
        for cath_line in cath_lines:
          cath_split = cath_line.split("\t")
          if cath_split[0] == target:
            # logger.info(cath_line)
            cath_id = cath_split[3]
            if cath_split[4].isdigit():
              start_uni = int(cath_split[4])
              if cath_split[5].strip("\n").isdigit():
                end_uni = int(cath_split[5].strip("\n"))

                for chain in uni_pdb:
                  # if entirely contained in pdb
                  if (start_uni in uni_pdb[chain] and 
                      end_uni in uni_pdb[chain]):
                    start_pdb = uni_pdb[chain][start_uni]
                    end_pdb = uni_pdb[chain][end_uni]
                    cath_numb[cath_id][chain] = (start_pdb, end_pdb)

                  # partial, no end
                  elif (start_uni in uni_pdb[chain] and 
                        (not end_uni in uni_pdb[chain])):
                    start_pdb = uni_pdb[chain][start_uni]
                    end_pdb = max(uni_pdb[chain].values())
                    cath_numb[cath_id][chain] = (start_pdb, end_pdb)
                  # partial, no start
                  elif ((not start_uni in uni_pdb[chain]) and 
                        (end_uni in uni_pdb[chain])):
                    start_pdb = min(uni_pdb[chain].values())
                    end_pdb = uni_pdb[chain][end_uni]
                    cath_numb[cath_id][chain] = (start_pdb, end_pdb)
              
              else:
                logger.info(str(cath_split[5]) + '!!!!')
            else:
                logger.info(str(cath_split[4]) + '!!!!')


        # logger.info(cath_numb)

        # sum of dictionaries
        arch_numb = dict(pfam_numb.items() + cath_numb.items())
        # logger.info('the sum of dics is' + str(arch_numb))

        # check for each residue if they are in the ranges?
        for het in het_ch_res:
          # loop over chain
          for chain in het_ch_res[het]:
            # each res we have in the list of interacting residues
            for res in het_ch_res[het][chain]:
              # loop over cath/pfam
              for domain in arch_numb:
                for chain in arch_numb[domain]:

                  if (int(res) < arch_numb[domain][chain][1] and 
                      int(res) > arch_numb[domain][chain][0]):
                    good_dom.append(domain)
              

          ######################################

      # still in the target loop, (for each drug) get the good targets
      # could be different for each target!

      # is there at least one good dom
      if good_dom:
        # get rid of duplicates
        good_dom = list(set(good_dom))
        # logger.info(' The good dom for target '+ str(target)
         # + ' are '+str(good_dom))
        # check each domain in our mapping file
        for arch in drug_arch_target[drug][target]:

          # if the arch is one of the good ones
          if arch in good_dom:
            drug_filt[drug][target][arch] = (
                                        drug_arch_target[drug][target][arch])

    # logger.info(pdb_uni)

  # filtered dictionaries with only entries that have domain interacting
  # with the drug
  return drug_filt

############################################################################




############################################################################
### STRUCT_MAPS
############################################################################
#take chembl or drugbank cluster (drg to het)
# 
# return dictionary drug:target:pdb:het
# struct_maps(
  # chembl_repo_map, chembl_cluster,
#                 uniprot_pdb_w_lig, pdb_cc_dic))

def struct_maps(repox, clust_het, uni_pdb, pdb_cc):

  # filtered version of big map, with only targets with structural info
  # drug: target: arch : target
  map1 = AutoVivification()

  # drug: target:pdb:het
  map2 = AutoVivification()

  # logger.info(uni_pdb)
  for drug in clust_het:
    # logger.info(drug)
    # dict for uniprot to pdbs
    lucky_uniprot = {}


    for protein in repox[drug]:
      # logger.info(protein)
      if protein in uni_pdb:
            # each pdb associated with them
            for pdb in uni_pdb[protein]:

              good_hets = []

              for listed_het in pdb_cc[pdb]:
                # logger.info(listed_het)
                if listed_het in clust_het[drug]:

                  good_hets.append(listed_het)


              # just in case
              if good_hets:

                good_hets = list(set(good_hets))
                # logger.info(pdb)
                # logger.info(good_hets)
                map1[drug][protein] = repox[drug][protein]
                
                map2[drug][protein][pdb] = good_hets



  return map1, map2

############################################################################



############################################################################
### PERCENT_IDENTITY
############################################################################
# from map get average percent identity after alignment
# map drug:target:arc:target
# return drug target: average identity score

def percent_identity(drug_targ_map):
  from Bio import pairwise2
  from Bio.SubsMat import MatrixInfo as matlist

  matrix = matlist.blosum62
  gap_open = -10
  gap_extend = -0.5

  # number_alignments = 0
  #overwrite map
  # drug_targ_map = {'CHEMBL98': {'Q92769': {'3.40.800.20': 
  # ['G4LW33', 'G4LW32','G4VAN3'], 'PF00850': ['G4LW33']}}}
 
  # dictionary for each target to store average identity score
  
  # targets_average_scores = {'P14324': [25, 28, 23, 15, 32],
   # 'P00519': [11, 24, 10, 25, 18, 20, 19, 21, 18, 24, 22, 23, 14, 
   # 13, 7, 18, 17, 6, 23, 23, 23, 16, 13, 16, 18, 20, 18, 14, 17, 15, 
   # 24, 25, 16, 18, 23, 20, 25, 21, 17, 24, 6, 18, 23, 25, 20, 23, 22,
    # 18, 25, 23, 24, 23, 15, 17, 18, 26, 23, 21, 23, 20, 25, 20, 21, 19, 
    # 3, 3, 22, 16, 23, 23, 17, 19, 24, 25, 24, 11, 23, 23, 21, 22, 16, 
    # 17, 15, 14, 21, 24, 28, 24, 25, 23, 14, 19, 22, 21, 22, 22, 21, 16, 
    # 20, 22, 19, 20, 16, 23, 20, 24, 23, 23, 23, 24, 23, 25, 14, 20, 20, 
    # 20, 5, 13, 19, 20, 23, 19, 19, 19, 13, 24, 19, 24, 25, 17, 22, 18, 
    # 20, 4, 11], 'P10827': [16, 15, 18, 24, 16, 15, 24, 23, 18, 14, 22, 
    # 17, 24, 24, 15, 25, 22, 23, 22, 22, 25], 'O76074': [21, 13, 23, 
    # 23, 26, 24, 26], 'P12931': [18, 24, 17, 20, 25, 24, 24, 26, 24, 
    # 24, 25, 23, 21, 19, 12, 25, 23, 10, 24, 21, 23, 23, 20, 23, 25, 
    # 23, 23, 21, 24, 25, 21, 23, 26, 25, 24, 24, 25, 24, 23, 22, 11, 
    # 24, 24, 21, 24, 23, 27, 24, 21, 22, 22, 24, 22, 26, 25, 20, 23, 
    # 26, 23, 25, 20, 24, 24, 27, 7, 7, 25, 22, 24, 19, 24, 26, 23, 
    # 21, 20, 18, 23, 17, 28, 23, 23, 23, 22, 22, 25, 18, 20, 25, 
    # 25, 24, 23, 24, 24, 23, 24, 25, 25, 22, 24, 27, 24, 26, 
    # 22, 21, 23, 22, 22, 17, 24, 20, 21, 22, 22, 24, 23, 25, 9, 
    # 21, 25, 24, 21, 23, 23, 24, 20, 19, 25, 23, 22, 23, 25, 24, 
    # 14, 8, 18], 'P10275': [21, 20, 22, 24, 22, 19, 17, 24, 22, 
    # 20, 22, 20, 24, 19, 19, 23, 16, 22, 25, 22, 23], 'P27487': 
    # [20, 14, 25, 23, 18, 17, 18, 18, 24, 24, 25, 26, 19, 25, 19, 
    # 26, 23, 20, 28, 22, 22, 17, 25, 16, 21, 16, 22, 24, 21], 
    # 'P03372': [18, 17, 20, 27, 19, 16, 23, 27, 20, 18, 24, 18, 
    # 26, 22, 17, 24, 21, 23, 25, 24, 24], 'P35968': [10, 23, 9, 
    # 25, 16, 18, 17, 20, 18, 22, 21, 21, 13, 11, 6, 17, 16, 5, 
    # 21, 22, 21, 17, 12, 15, 17, 19, 17, 13, 15, 14, 25, 23, 15, 
    # 16, 21, 18, 23, 19, 17, 24, 5, 16, 22, 24, 18, 22, 21, 16, 
    # 24, 23, 22, 23, 14, 16, 16, 24, 22, 20, 22, 18, 25, 18, 19, 
    # 19, 3, 3, 24, 14, 23, 24, 17, 18, 22, 24, 25, 10, 22, 24, 
    # 21, 22, 15, 17, 14, 13, 20, 24, 24, 22, 23, 22, 13, 17, 22, 
    # 19, 20, 22, 20, 14, 18, 21, 18, 19, 14, 24, 19, 23, 25, 
    # 24, 21, 24, 25, 26, 13, 18, 19, 19, 4, 13, 18, 18, 23, 18, 
    # 18, 18, 12, 25, 17, 21, 25, 15, 21, 17, 21, 3, 10]}

  for drug in drug_targ_map:
  # for drug in ['CHEMBL98', 'CHEMBL973']:
    logger.info(drug)
    for targ in drug_targ_map[drug]:
      # check it in not in there already
      
      target_identity_scores = []

      if targ not in targets_average_scores:
        logger.info(targ)
        # list identity score for each target
        # logger.info(targ)

        # write fasta drug target
        try:
          # get handle
          handle = ExPASy.get_sprot_raw(targ)
          # swissprot read
          targ_record = SwissProt.read(handle)

          # sequence
          targ_seq = targ_record.sequence
          
        except HTTPError:
          logger.error('We cannot access ExPASy!')
          logger.warning('We are exiting the program.')
          sys.exit()


        # logger.info(targ_seq)
        
        # get uniwue list of proteins to align to this target
        # because there will be duplicates!
        protein_list = []
        
        for arch in drug_targ_map[drug][targ]:
          for schi in drug_targ_map[drug][targ][arch]:
            protein_list.append(schi)


        logger.info(len(protein_list))
        protein_list = list(set(protein_list))
        logger.info(len(protein_list))  

        # each schisto protein
        for prot in protein_list:

          try:
              # get handle
              handle = ExPASy.get_sprot_raw(prot)
              # swissprot read
              record = SwissProt.read(handle)

              # sequence
              seq = record.sequence
              

          except HTTPError:
            logger.error('We cannot access ExPASy!')
            logger.warning('We are exiting the program.')
            sys.exit()

          # logger.info(seq)
          
          # alns = pairwise2.align.globalds(targ_seq, seq, matrix, 
            # gap_open, gap_extend)
          alns = pairwise2.align.globalxx(targ_seq, seq)

          percent = int(alns[0][2]) * 100 / int(alns[0][4])
          # logger.info(percent)
          target_identity_scores.append(percent)

      # if target is already in there, do nothing
      else:
        pass
      # logger.info(target_identity_scores)
      # target_average = sum(target_identity_scores)/
      # len(target_identity_scores)
      # logger.info('the average for one target is'+str(target_average))
      # logger.info(target_average)

      # only if there is something in the list
      if target_identity_scores:
        # add list to the dictionary
        targets_average_scores[targ] = target_identity_scores
        logger.info(targets_average_scores)

  
  return targets_average_scores



############################################################################




############################################################################
### MAIN FUNCTION
############################################################################

def main():
  # time
  start_time = datetime.now()

  # pipeline step counter
  step = 0
  
  # header
  logger.info('********************** DRUG_REPO.PY ' +
              '**********************')
  # greeting
  logger.info("Hi " + c.your_name + 
            ", you are running drug_repo.py for drug repositioning.")


  # check number of steps
  if c.steps not in range(0,11):
    logger.error('You must select a step number between 0 and 10!')
    logger.info('Please check config.py')
    logger.warning('The program is aborted.')
    sys.exit()
  else:
    pass

 
  # check species list
  if not c.taxa:
    logger.error('The list of species is empty!')
    logger.info('Please check config.py')
    logger.warning('The program is aborted.')
    sys.exit()
  else:
    # create list of species name
    species_list = taxa_to_species(c.taxa, c.spec_list)

    # create string of species
    species_string = lst_to_string(species_list)

    # check the species string is not empty
    if species_string == '':
      logger.error('We could not find the species you have selected in the' +
                    ' mapping file!')
      logger.info('Please check config.py')
      logger.warning('The program is aborted.')
      sys.exit()
    else:
      pass

  # check datasets (chembl vs drugbank)
  dataset_list = []
  for dataset in c.sets:
    if dataset in c.dataset_dic:
      # add dataset name to list
      dataset_list.append(c.dataset_dic[dataset])

    else:
      logger.error('The list of datasets contains unallowed letters!')
      logger.info('Please check config.py')
      logger.warning('The program is aborted.')
      sys.exit()

  if dataset_list:
    dataset_string = lst_to_string(dataset_list)

  else:
    logger.error('The list of datasets is empty!')
    logger.info('Please check config.py')
    logger.warning('The program is aborted.')
    sys.exit()



  logger.info("You have selected to run the pipeline up to step number " + 
              str(c.steps) + ", against dataset(s) from " +
              dataset_string + ", and to investigate species " + 
              species_string + ".")

  logger.info("If that does not sound right" +
              ", please customise your settings in config.py " +
              "before proceeding. Refer to README for more info. " +
              "If you are all set, let's do some mapping!")
  


  #################################
  ### STEP 1: FIND DRUG TARGETS ###
  #################################
  step = step + 1

  if c.steps > (step-1):

    logger.info('------------------------- STEP ' + str(step) + ' ' +
                '-------------------------')
    logger.info('We wish to take the ' + dataset_string + 
                ' dataset(s), map the drug ids to their target ids ' +
                'and obtain a list of unique drug targets.')
    
    # SET A (chembl)
    if 'A' in c.sets:
      logger.info('We are processing the ' + str(c.dataset_dic['A']) + 
                  'input file, ' + c.chembl_input)

      # generate chembl dictionary
      chembl_dic = run_or_pickle("1_chembl_dic", process_chembl, 
                                  c.chembl_input)

      # logger.info(chembl_dic)

      # get list of uniprot ids from chembl_dic
      chembl_uniprot_list = run_or_pickle("1_chembl_uniprot_list",
                                          flatten_dic, chembl_dic, "values")

      logger.info(str(len(chembl_dic)) + ' ChEMBL drugs could be mapped to ' +
                   str(len(chembl_uniprot_list)) + ' unique UniProt ID.')
      
      taxa_chembl = run_or_pickle("1_human_chembl", expasy_dic, 
                                    chembl_uniprot_list, "taxa")
      human_chembl = taxa_chembl['9606']

      percent_chembl_human = round(float(len(human_chembl)) / 
                            float(len(chembl_uniprot_list)) * 100)

      logger.info('Of these targets, ' + str(len(human_chembl)) + 
                  ' are human proteins (the ' + str(percent_chembl_human) + 
                  ' %).')

    else:
      chembl_uniprot_list = []


    # if set B (drugbank)
    if 'B' in c.sets:
      logger.info('We are processing the ' + str(c.dataset_dic['B']) + 
                  'input file, ' + c.drugbank_input)

      # generate drugbank_dictionary
      drugbank_dic = run_or_pickle("1_drugbank_dic", process_drugbank, 
                      c.drugbank_input)

      # logger.info(drugbank_dic)
      #logger.info(len(drugbank_dic))

      # get list of uniprot ids from drugbank
      drugbank_uniprot_list = run_or_pickle("1_drugbank_uniprot_list",
                                        flatten_dic, drugbank_dic, "values")
      

      logger.info(str(len(drugbank_dic)) + 
                  ' DrugBank drugs could be mapped to ' +
                  str(len(drugbank_uniprot_list)) + ' unique UniProt ID.')

    else:
      drugbank_uniprot_list = []


    # merge lists and rm duplicates
    uniprot_list = run_or_pickle("1_uniprot_list", merge_lists,
                                chembl_uniprot_list, drugbank_uniprot_list)

    logger.info('We have obtained a unique list of ' + 
              str(len(uniprot_list)) + ' drug targets.')

    #logger.info(uniprot_list)
    #uniprot_list = ['Q09428', 'Q5BVT7', 'P00519']
    # reviewed_targets = run_or_pickle("1_reviewed_targets", expasy_filter, 
    #                                   uniprot_list, "reviewed") 



    # dictionary of taxa codes vs lists of uniprot
    # to evaluate how many human targets and other species
    # taxa_targets = run_or_pickle("1_human_targets", expasy_dic, 
                                    # uniprot_list, "taxa")


    # percent_human = round(float(len(taxa_targets['9606'])) / 
    #                       float(len(uniprot_list)) * 100)
    
    # ### OVERWRITE UNIPROT_LIST WITH MADE-UP LIST
    # #uniprot_list = ['Q4JEY0', 'P68363', 'P10613','P18825']
    # ###


    # logger.info('Of these targets, ' + str(len(taxa_targets['9606'])) + 
                # ' are human proteins (the ' + str(percent_human) + ' %).')
    logger.info('----------------------------------------------------------')
    

  else:

    logger.info('We have skipped STEP ' + str(step))  


  ############################
  ### STEP 2: FIND ARCH    ###
  ############################
  step = step + 1

  if c.steps > (step-1):
    logger.info('------------------------- STEP ' + str(step) + ' ' +
                '-------------------------')

    logger.info('We wish to build dictionaries of ' +
                'the targets\' Uniprot ids to CATH/Uniprot ids.')

    # splitting chembl_het_map in chunks
    uniprot1,uniprot2,uniprot3, uniprot4, uniprot5 = \
                        zip(*izip_longest(*[iter(uniprot_list)]*5))

    logger.info('We have split the uniprot list in chuncks of ' +
                str(len(uniprot1)) + ', ' + str(len(uniprot2)) + ', ' + 
                str(len(uniprot3)) + ', ' + str(len(uniprot4)) + ', ' +
                str(len(uniprot5)) )



    # run or pickle uniprot_to_arch to retrieve cath domain architectures
    cath_dic1 = run_or_pickle("2_cath_dic1", uniprot_to_arch, uniprot1, 
                            "cath")

    cath_dic2 = run_or_pickle("2_cath_dic2", uniprot_to_arch, uniprot2, 
                            "cath")
    cath_dic3 = run_or_pickle("2_cath3_dic3", uniprot_to_arch, uniprot3, 
                            "cath")
    cath_dic4 = run_or_pickle("2_cath3_dic4", uniprot_to_arch, uniprot4, 
                            "cath")
    cath_dic5 = run_or_pickle("2_cath3_dic5", uniprot_to_arch, uniprot5, 
                            "cath")


    cath_dic = dict(cath_dic1.items() + 
                  cath_dic2.items() + cath_dic3.items() + cath_dic4.items() +
                  cath_dic5.items())


    #logger.info(cath_dic)
    # generate list, flatten it and rm duplicates
    cath_list = run_or_pickle("2_cath_list", flatten_dic, cath_dic, "values")



    logger.info('We have mapped ' + str(len(cath_dic)) + ' uniprot ids to ' +
                str(len(cath_list)) + ' CATH ids.')

    # run or pickle uniprot_to_arch to retrieve pfam domain architectures
    pfam_dic1 = run_or_pickle("2_pfam_dic1", uniprot_to_arch, uniprot1, 
                            "pfam")

    pfam_dic2 = run_or_pickle("2_pfam_dic2", uniprot_to_arch, uniprot2, 
                            "pfam")

    pfam_dic3 = run_or_pickle("2_pfam_dic3", uniprot_to_arch, uniprot3, 
                            "pfam")
    pfam_dic4 = run_or_pickle("2_pfam_dic4", uniprot_to_arch, uniprot4, 
                            "pfam")

    pfam_dic5 = run_or_pickle("2_pfam_dic5", uniprot_to_arch, uniprot5, 
                            "pfam")



    pfam_dic = dict(pfam_dic1.items() + 
                  pfam_dic2.items() + pfam_dic3.items() + 
                  pfam_dic4.items() + pfam_dic5.items())


    #logger.debug(len(pfam_dic))

    # generate list, flatten it and rm duplicates
    pfam_list = run_or_pickle("2_pfam_list", flatten_dic, pfam_dic, "values")


    logger.info('We have mapped ' + str(len(pfam_dic)) + ' uniprot ids to ' +
                str(len(pfam_list)) + ' Pfam ids.')

    logger.info('----------------------------------------------------------')


  else:
    logger.info('We have skipped STEP ' + str(step))

  
  ####################################
  ### STEP 3: FIND SCHISTO TARGETS ###
  ####################################
  step = step + 1

  if c.steps > (step-1):
    logger.info('------------------------- STEP ' + str(step) + ' ' +
                '-------------------------')

    logger.info('We wish to map the CATH/Pfam ids ' +
                'to UniProt ids of the species ' + species_string + '.')
    
   
    if len(cath_list) < 1000:
      # call archindex on cath values to find the ones from schisto
      uniprot_schisto_cath_dic = run_or_pickle("3_uniprot_schisto_cath_dic", 
                                            arch_to_uniprot, cath_list, 
                                            "cath")

    else:

      # splitting cath_list in chunks
      cath_list1, cath_list2, cath_list3 = \
                        zip(*izip_longest(*[iter(cath_list)]*3))
      # logger.info(len(cath_list2))

      uniprot_schisto_cath_dic1 = run_or_pickle("3_uniprot_schisto_cath_dic1", 
                                            arch_to_uniprot, cath_list1, 
                                            "cath")
      uniprot_schisto_cath_dic2 = run_or_pickle("3_uniprot_schisto_cath_dic2", 
                                            arch_to_uniprot, cath_list2, 
                                            "cath")
      uniprot_schisto_cath_dic3 = run_or_pickle("3_uniprot_schisto_cath_dic3", 
                                            arch_to_uniprot, cath_list3, 
                                            "cath")

      uniprot_schisto_cath_dic = dict(uniprot_schisto_cath_dic1.items() + 
                  uniprot_schisto_cath_dic2.items() + 
                  uniprot_schisto_cath_dic3.items())


    # generate list, flatten it and rm duplicates
    uniprot_schisto_cath_list = run_or_pickle("3_uniprot_schisto_cath_list",
                                flatten_dic, uniprot_schisto_cath_dic, 
                                "values")

    logger.info('We have mapped ' + str(len(uniprot_schisto_cath_dic)) + 
                ' CATH ids to ' + str(len(uniprot_schisto_cath_list)) +
                ' Uniprot ids.')

    # call archindex on pfam values to find ones from schisto
    uniprot_schisto_pfam_dic = run_or_pickle("3_uniprot_schisto_pfam_dic", 
                                            arch_to_uniprot, pfam_list, 
                                            "pfam")
    #logger.debug(len(uniprot_schisto_pfam_dic))



    # generate list, flatten it and rm duplicates
    uniprot_schisto_pfam_list = run_or_pickle("3_uniprot_schisto_pfam_list",
                                flatten_dic, uniprot_schisto_pfam_dic, 
                                "values")

    # logger.info(uniprot_schisto_pfam_dic)
    logger.info('We have mapped ' + str(len(uniprot_schisto_pfam_dic)) + 
                ' pfam ids to ' + str(len(uniprot_schisto_pfam_list)) +
                ' Uniprot ids.')

    # merge and rm duplicates
    # this is total list of unique schisto uniprot ids
    uniprot_schisto_list = run_or_pickle("3_uniprot_schisto_list", 
                          merge_lists, uniprot_schisto_cath_list, 
                          uniprot_schisto_pfam_list)
    
    logger.info('In total, we have identified ' + 
                str(len(uniprot_schisto_list)) + 
                ' unique targets that point to known drugs.')


    # filter list for only reviewed ones
    # uniprot_schisto_filt = run_or_pickle("3_uniprot_schisto_filt",
    #                                       expasy_filter,
    #                                       uniprot_schisto_list, "reviewed")

    # ### OVERWRITE FILTERED LIST
    # #uniprot_schisto_filt = ['P33676']
    # ###

    # logger.info('Of those targets, the reviewed Uniprot entries are ' + 
    #             str(len(uniprot_schisto_filt)) +  '.')
    
    logger.info('----------------------------------------------------------')

  else:
    logger.info('We have skipped STEP ' + str(step))


  ############################
  ### STEP 4 GENERATE MAPS ###
  ############################
  step = step + 1

  if c.steps > (step-1):
    logger.info('------------------------- STEP ' + str(step) + ' ' +
                '-------------------------')

    logger.info('We wish to create dictionaries that collect all the' +
                ' mapping so far.')


    # SET A (chembl)
    if 'A' in c.sets:

      # generate big map for chembl drugs
      # drug: target: arch: targ
      chembl_repo_map = run_or_pickle("4_chembl_repo_map", chembl_repo, 
                                      chembl_dic, cath_dic,
                                      uniprot_schisto_cath_dic, pfam_dic, 
                                      uniprot_schisto_pfam_dic)
      # logger.debug(len(chembl_repo_map))
      chembl_repo_schisto_list = flatten_dic(chembl_repo_map, 'values_4')

      # drugtarg = []

      # for drug in chembl_repo_map:
      #   for tar in chembl_repo_map[drug]:
      #     drugtarg.append(tar)

      # drugtarg = list(set(drugtarg))
      # logger.info(len(drugtarg))


      logger.info('We have built the ChEMBL map, mapping ' +
                  str(len(chembl_repo_map)) + ' ChEMBL drugs to ' +
                  str(len(chembl_repo_schisto_list)) + ' unique ' +
                  species_string + ' targets.')

      # list of drugs that are in the map, to be used in part 6
      #chembl_repo_drug_list = chembl_repo_map.keys()
      chembl_repo_drug_list = flatten_dic(chembl_repo_map, 'keys')


      # obtain filtered mapping dictionary for filtered entries
      # chembl_schisto_filt_map = run_or_pickle("4_chembl_schisto_filt_map",
      #                                       filt_schisto_map, chembl_repo_map,
      #                                       uniprot_schisto_filt)

    else:
      pass


    # SET B (drugbank)
    if 'B' in c.sets:

      # generate big map for drugbank drugs
      drugbank_repo_map = run_or_pickle("4_drugbank_repo_map", drugbank_repo, 
                                        drugbank_dic, cath_dic,
                                        uniprot_schisto_cath_dic, pfam_dic, 
                                        uniprot_schisto_pfam_dic)
      


      logger.info('We have built the DrugBank map, mapping ' +
                  str(len(drugbank_repo_map)) + 
                  ' DrugBank drugs to potential ' +
                  species_string + ' targets.')
      # logger.info(drugbank_repo_map)

      # list of drugs that are in the map, to be used in part 6
      # below old one, had white spaces!
      #drugbank_repo_drug_list = drugbank_repo_map.keys()
      #logger.debug(drugbank_repo_drug_list)
      # new list, no white spaces!
      drugbank_repo_drug_list = flatten_dic(drugbank_repo_map, 'keys')

      # filtered ap for reviewed entries!

      # drugbank_schisto_filt_map = run_or_pickle("4_drugbank_schisto_filt_map",
      #                                           filt_schisto_map, 
      #                                           drugbank_repo_map,
      #                                           uniprot_schisto_filt)

    else:
      pass

    logger.info('----------------------------------------------------------')


  else:
    logger.info('We have skipped STEP ' + str(step))


  #################################################
  ### STEP 5 PDB TO HET GROUPS                  ###
  #################################################
  step = step + 1

  if c.steps > (step-1):
    logger.info('------------------------- STEP ' + str(step) + ' ' +
                '-------------------------')

    logger.info('We wish to map all available pdb structures ' +
                'to the Het groups the contain, and then filter out ' +
                'the Het groups contained in ' + c.pointless_het + 
                ', a list of ions, metals, peptidic ligands, etc.')
    # make dictionary of pdb to ligands
    pdb_lig_dic = run_or_pickle("5_pdb_lig_dic", lst_dic, c.pdb_lig)

    logger.info('We made a dictionary of '+ str(len(pdb_lig_dic)) + 
              ' pdb entries mapped to their ligand identifiers.')

    # make list of ccs to ignore
    pointless_het = run_or_pickle("5_pointless_het", csv_to_lst,
                                      c.pointless_het)
    logger.info("The list of ligands we wish to ignore " +
                "contains " + str(len(pointless_het)) + " ligands.")


    # filter dictionary pdb to lig, excluding lig that are in the 'pointless'
    # list
    pdb_lig_pointless_dic = run_or_pickle("5_pdb_lig_pointless_dic",
                                      exclude_values_from_dic, pdb_lig_dic, 
                                      pointless_het, "exclude")
    #logger.debug(len(pdb_lig_pointless_dic))
    
    # regular expression for string containing at least one dash
    contains_dash = re.compile('.*-.*')

    # second filter, to eliminate those with dash
    # this is filtered dic of all useful pdbs (with useful ligands!)
    pdb_lig_filt_dic = run_or_pickle("5_pdb_lig_filt_dic",
                                      exclude_values_from_dic, 
                                      pdb_lig_pointless_dic, 
                                      contains_dash, "nomatch")
    #logger.info(pdb_lig_filt_dic)
    
    # list of 'acceptable' pdbs (with useful ligands) from dic
    pdb_w_lig_list = run_or_pickle("5_pdb_w_lig_list", flatten_dic, 
                                  pdb_lig_filt_dic, "keys")
    # the length here is obviously the same as the length of dic!
    #logger.info(pdb_w_lig_list)

    # list of 'acceptable' ligands from dic
    filtered_ligs = run_or_pickle("5_filtered_ligs", flatten_dic,
                                  pdb_lig_filt_dic, "values")
    #logger.debug(filtered_ligs)


    logger.info('We have excluded the pdb entries that only have ' +
                'ligands from such list, ' + 
                'to obtain ' + str(len(pdb_lig_filt_dic)) +
                ' pdb entries, mapped to a total of ' +
                str(len(filtered_ligs)) + ' unique ligands.')

    logger.info('----------------------------------------------------------')

  else:
    logger.info('We have skipped STEP ' + str(step))

  ##################################################
  ### STEP 6 DRUG TARGETS WITH STRUCTURAL INFO   ###
  ##################################################
  step = step + 1

  if c.steps > (step-1):
    logger.info('------------------------- STEP ' + str(step) + ' ' +
                '-------------------------')

    logger.info('We wish to collect all the drug targets that ' +
                'point to some repositioning target, point them to ' +
                'the available pdb structures (using ' + c.uniprot_pdb + 
                  '), filter them according to the map obtained in Part 4 ' +
                  'and extract the ligands (using ' + c.cc_smi + ').')
    # obtain list of targets from drugbank_repo_map and chembl_repo_map
    # these are all the uniprot values that are targets of our potential
    # drug repo candidates


    # SET A (chembl)
    if 'A' in c.sets:
      chembl_drug_targ = run_or_pickle("6_chembl_drug_targ", 
                                    list_second_level_dic, chembl_repo_map)
    else:
      chembl_drug_targ = []

    # SET B (drugbank)
    if 'B' in c.sets:
      drugbank_drug_targ = run_or_pickle("6_drugbank_drug_targ", 
                                    list_second_level_dic, drugbank_repo_map)
    else:
      drugbank_drug_targ = []

    

    tot_drug_targ = run_or_pickle("6_tot_drug_targ", merge_lists, 
                                  chembl_drug_targ, drugbank_drug_targ)

    logger.info('Overall, we have ' + str(len(tot_drug_targ)) + 
                ' drug targets that could be mapped to some schisto target.')

    # make dictionary uniprot to pdb
    uniprot_pdb_dic = run_or_pickle("6_uniprot_pdb_dic", csv_to_dic, 
                                    c.uniprot_pdb)
    #logger.debug(uniprot_pdb_dic)

    # this is dictionary of drug targets that have at least one pdb structure
    uniprot_filt = run_or_pickle("6_uniprot_filt", filter_dic_from_list, 
                                uniprot_pdb_dic, tot_drug_targ)



    logger.info('Of those targets, ' + str(len(uniprot_filt)) + 
                ' have at least one pdb structure associated to them.')
    # apply filter
    #logger.debug(uniprot_filt)
   
    #####
    # this would be alternative method for finding entries with pdb
    # BUT! it returns larger list (eg 2751 instead of 2711) 
    # because it includes
    # pdb ids that point to model structures, not accepted in the pdb anymore
    # # obtain filtered list of drug targets that have associated pdbs
    # uniprot_filt_pdb = run_or_pickle("5_uniprot_filt_pdb", expasy_filter, 
    #                                 uniprot_list, "pdb")
    #####


    # take dic of drug target uniprot to pdb and keep only ones that 
    # are in in the 'acceptable' pdb list
    uniprot_pdb_w_lig = run_or_pickle("6_uniprot_pdb_w_lig", 
                                      exclude_values_from_dic, uniprot_filt,
                                      pdb_w_lig_list, "include")
    #logger.info(uniprot_pdb_w_lig)


    ###
    # get list of uniprot from dic above
    uniprot_w_lig_list = uniprot_pdb_w_lig.keys()
    #logger.info(len(uniprot_w_lig_list))

    ###

    # get the pdb list from the dic above
    pdb_w_lig = run_or_pickle("6_pdb_w_lig", list_second_level_dic,
                                  uniprot_pdb_w_lig)


    
    logger.info('Of those, ' + str(len(uniprot_pdb_w_lig)) +
                ' have at least one pdb structure ' +
                'in complex with a small molecule associated to them. ' +
                'The unique pdb structures are ' + str(len(pdb_w_lig)) + '.')


    # now filter the pdb_lig_filt_dic, to obtain the pdbs we want
    # { PDBID: [list of CC]}
    # this will be the dic we refer to later!
    pdb_cc_dic = run_or_pickle("6_pdb_cc_dic", filter_dic_from_list, 
                                  pdb_lig_filt_dic, pdb_w_lig)


    # finally obtain the list of cc we need! - the ones that are in the pdbs
    # this is the list of cc we need to try and match to the drugs cc!!
    # will be around 8435 in the list
    cc_list = run_or_pickle("6_cc_list", list_second_level_dic,
                            pdb_cc_dic)
    logger.info('We have extracted the chemical components from the pdbs, ' +
                'for a total of ' + str(len(cc_list)) + 
                ' chemical components')

    # get cc to smiles dictionary
    cc_smiles = run_or_pickle("6_cc_smiles", smi_to_dic, c.cc_smi, 1, 0)

    #logger.info(len(cc_smiles))


    # dic of cc we ar interested in, mapped to their smiles
    cc_smi_filt = run_or_pickle("6_cc_smi_filt", filter_dic_from_list, 
                                cc_smiles, cc_list)
    
    logger.info('We have mapped ' + str(len(cc_smi_filt)) + 
                ' of these chemical components to their smiles.')
    
    logger.info('----------------------------------------------------------')

  else:
    logger.info('We have skipped STEP ' + str(step))

  ####################################
  ### STEP 7 COMPOUND CLUSTERING   ###
  ####################################
  step = step + 1

  if c.steps > (step-1):
    logger.info('------------------------- STEP ' + str(step) + ' ' +
                '-------------------------')

    logger.info('We wish to take the drugs ' +
                'from the mapping and cluster them against ' +
                'the chemical components extracted from the pdb structures.')

    # SET A (chembl)
    if 'A' in c.sets:

      # total chembl drugs to smiles dictionary - 10406 chembl drugs
      chembl_id_smi_dic = run_or_pickle("7_chembl_id_smi_dic", txt_to_dic, 
                                        c.chembl_input, "CHEMBL_ID",
                                        "CANONICAL_SMILES")
      #logger.debug(len(chembl_id_smi_dic))
      

      # filter dictionary to only drugs that are in chembl_repo_drug_list
      # these are all chembl drugs (783) that are in the map
      chembl_id_smi_filt = run_or_pickle("7_chembl_id_smi_filt", 
                                          filter_dic_from_list, 
                                          chembl_id_smi_dic,
                                          chembl_repo_drug_list)
      
      logger.info('We have mapped the ' + str(len(chembl_id_smi_filt)) +
                  ' ChEMBL drugs to their smiles.')


      # filter chembl_dic to only the 783 drugs, using chembl_repo_drug_list
      chembl_dic_mapped_drugs = filter_dic_from_list(chembl_dic, 
                                chembl_repo_drug_list)
      #logger.info(len(chembl_dic_mapped_drugs))
      # filter out the uniprots, using uniprot_w_lig_list
      chembl_dic_uni_drugs = exclude_values_from_dic(chembl_dic_mapped_drugs, 
                            uniprot_w_lig_list, "include")
      #logger.info(len(chembl_dic_uni_drugs))
      chembl_uni_drugs_list = chembl_dic_uni_drugs.keys()
      #logger.info(len(chembl_uni_drugs_list))
      # filter chembl_id_smi_filt to what obtained above
      chembl_id_smi_opt = filter_dic_from_list(chembl_id_smi_filt, 
                          chembl_uni_drugs_list) 
      #logger.info(len(chembl_id_smi_opt))
      

      # obtain drug to cc dictionary, merging three dics
      # this is all the drugs in the map, pointing to the cc in the pdbs of
      # of their targets
      chembl_to_cc = merge_dic(chembl_dic,uniprot_filt, pdb_cc_dic)
      #logger.info(len(chembl_to_cc))


      logger.info('We have filtered out the ChEMBL drugs that ' +
                  'do not point to a crystal structure in complex with ' +
                  'a small molecule, to obtain ' + 
                    str(len(chembl_id_smi_opt)) +
                  ' ChEMBL drugs mapped to their smiles' +
                  '; these will be clustered.')


      #####################
      #openbabel conversion - not necessary for now
      #####################
      # # cc: create file with smiles to feed to openbabel
      # dic_to_txt(cc_smi_filt, '6_cc_smi_filt.smi')

      # # ligands converted to sdf
      # # babel_smi_to_sdf('test.smi','test.sdf')

      # # ligands converted to sdf 3d and no hydrogens
      # #babel_smi_to_sdf('test.smi','test_noh.sdf')

      # # convert smiles to sdf
      # babel_smi_to_sdf('6_cc_smi_filt.smi','6_cc_smi_filt.sdf')
      #######################

      

      # # CHEMBL CLUSTERING (takes an hour approx)
      # rn clustering with Tanimoto similarity threshold
      # thresholds 1, 0.9, 0.8, 0.7 are also written to output:

      chembl_cluster = run_or_pickle("7_chembl_cluster", run_smsd, 
                                    chembl_id_smi_opt, cc_smi_filt,
                                    "pair_2dic", c.sim_threshold, 
                                    chembl_to_cc)
      # logger.info(chembl_cluster)
      # move output file to current dir
      mv_file(c.smsd_path, 'smsd_run_pair_2dic.txt', 
              c.chembl_clust_sim_scores)

      logger.info('We have clustered the ChEMBL drugs, to obtain ' + 
                  str(len(chembl_cluster)) + ' drugs mapped to at least ' +
                  'a chemical component with Tanimoto similarity above ' +
                  str(c.sim_threshold) + 
                  ' (other similarity thresholds written to ' +
                    str(c.chembl_clust_sim_scores) + ').')

      # get the list of drugs from cluster dic
      chembl_cluster_list = flatten_dic(chembl_cluster, "keys")
      # logger.info(chembl_cluster_list)

      # write filtered txt csv file to be imported in excel
      filter_txt('chembl_drugs.txt', c.chembl_cluster, 'CHEMBL_ID', 
                 chembl_cluster_list)

      
      logger.info('We have written to file ' + str(c.chembl_cluster) +
                  ' the info from ChEMBL regarding the clustered drugs' +
                  ', to be imported in excel.')



    
      # map chembl drugs to target to pdb to het
      chembl_struct_map, chembl_het_map = (
                                  struct_maps(chembl_repo_map, chembl_cluster,
                                              uniprot_pdb_w_lig, pdb_cc_dic))
      
      # chembl_struct_map ---> drug:target:arch:schisto target
      # chembl_het_map ----> drug: target: pdb: het

      # get list of schisto targets
      schis_targ = []
      for drug in chembl_struct_map:
        for tar in chembl_struct_map[drug]:
          for ar in chembl_struct_map[drug][tar]:
            for sch in chembl_struct_map[drug][tar][ar]:
              schis_targ.append(sch)

      schis_targ = list(set(schis_targ))
      
      logger.info('At this stage we have obtained ' + 
                  str(len(chembl_struct_map)) + ' drugs, mapped to ' +
                  str(len(schis_targ)) +
                  ' ' + species_string + ' targets.')

      # at this stage it is all the ones with struct info, but all domains!!

    else:
      pass


    # SET B (drugbank)
    if 'B' in c.sets:

      # drugbank drugs to smiles dictionary (total 6799 drugs mapped to smiles)
      drugbank_id_smi_dic = run_or_pickle('8_drugbank_id_smi_dic', 
                                          sdf_to_dic, c.drugbank_sdf, 
                                          'DATABASE_ID', 'SMILES')
      # logger.info(drugbank_id_smi_dic)
      

      # filter dictionary to only drugs that in the drugbank_repo_drug_list
      drugbank_id_smi_filt = run_or_pickle("8_drugbank_id_smi_filt", 
                                          filter_dic_from_list, 
                                          drugbank_id_smi_dic,
                                          drugbank_repo_drug_list)
      
      #logger.info(len(drugbank_id_smi_filt))
      #logger.info(drugbank_repo_drug_list)

     
      logger.info('We have mapped ' + str(len(drugbank_id_smi_filt)) +
                  ' DrugBank drugs to their smiles.')

      # obtain drug to cc dictionary, merging three dics
      drugbank_to_cc = merge_dic(drugbank_dic,uniprot_filt, pdb_cc_dic)
      # logger.info(drugbank_to_cc)

      
      logger.info('We have filtered the drugs, to obtain ' + 
                  str(len(drugbank_to_cc)) +
                  ' DrugBank drugs that will be clustered.')


      # splitting drugbank_to_cc into 5 chunks
      items1,items2,items3, items4, items5 = \
                          zip(*izip_longest(*[iter(drugbank_to_cc.items())]*5))
      d1 = dict(item for item in items1 if item is not None)
      d2 = dict(item for item in items2 if item is not None)
      d3 = dict(item for item in items3 if item is not None)
      d4 = dict(item for item in items4 if item is not None)
      d5 = dict(item for item in items5 if item is not None)
      
      logger.info('We have split the drugbank entries into chunks of ' +
                  str(len(d1)) + ', ' + str(len(d2)) + ', ' +
                  str(len(d3)) + ', ' + str(len(d4)) + ' and ' + 
                  str(len(d5)) + ', for easier processing.')


      # 1st
      logger.info('We are processing the first chunk.')
      db1_cluster = run_or_pickle("8_db1_cluster", run_smsd, 
                                    drugbank_id_smi_filt, cc_smi_filt,
                                    "pair_2dic", c.sim_threshold, d1)
    
      mv_file(c.smsd_path, 'smsd_run_pair_2dic.txt', '8_db1_cluster.txt')
      #logger.info(len(drugbank_cluster))
      #2nd
      logger.info('We are processing the second chunk.')
      db2_cluster = run_or_pickle("8_db2_cluster", run_smsd, 
                                    drugbank_id_smi_filt, cc_smi_filt,
                                    "pair_2dic",c.sim_threshold , d2)
      
      mv_file(c.smsd_path, 'smsd_run_pair_2dic.txt', '8_db2_cluster.txt')

      #3rd
      logger.info('We are processing the third chunk.')
      db3_cluster = run_or_pickle("8_db3_cluster", run_smsd, 
                                    drugbank_id_smi_filt, cc_smi_filt,
                                    "pair_2dic", c.sim_threshold, d3)
      
      mv_file(c.smsd_path, 'smsd_run_pair_2dic.txt', '8_db3_cluster.txt')

      #4th
      logger.info('We are processing the fourth chunk.')
      db4_cluster = run_or_pickle("8_db4_cluster", run_smsd, 
                                    drugbank_id_smi_filt, cc_smi_filt,
                                    "pair_2dic", c.sim_threshold, d4)
      
      mv_file(c.smsd_path, 'smsd_run_pair_2dic.txt', '8_db4_cluster.txt')

      #5th
      logger.info('We are processing the fifth chunk.')
      db5_cluster = run_or_pickle("8_db5_cluster", run_smsd, 
                                    drugbank_id_smi_filt, cc_smi_filt,
                                    "pair_2dic", c.sim_threshold, d5)
      
      mv_file(c.smsd_path, 'smsd_run_pair_2dic.txt', '8_db5_cluster.txt')

      # sum of all the 5 dics!
      drugbank_cluster = dict(db1_cluster.items() + db2_cluster.items() + 
                    db3_cluster.items() + db4_cluster.items() + 
                    db5_cluster.items())

      logger.info('We have clustered the DrugBank drugs, to obtain ' + 
                  str(len(drugbank_cluster)) + ' drugs mapped to at least ' +
                  'a chemical component with Tanimoto similarity above ' +
                  str(c.sim_threshold) + 
                  ' (other similarity thresholds written to file).')
    else:
      pass



    logger.info('----------------------------------------------------------')

  else:
    logger.info('We have skipped STEP ' + str(step))

  ####################################
  ### STEP 8 DOMAIN-RES FILTER   ###
  ####################################

  step = step + 1

  if c.steps > (step-1):
    logger.info('------------------------- STEP ' + str(step) + ' ' +
                '-------------------------')
    
    logger.info('We are now filtering the results to obtain ' +
                'a refined dictionary of drugs mapped to potential ' + 
                species_string + ' targets. We will only include the ' +
                'domains that are known to interact with drug' +
                ', or a close analogue.')
    


    # splitting chembl_het_map in chunks
    items1,items2,items3 = \
                        zip(*izip_longest(*[iter(chembl_het_map.items())]*3))
    ch_het1 = dict(item for item in items1 if item is not None)
    ch_het2 = dict(item for item in items2 if item is not None)
    ch_het3 = dict(item for item in items3 if item is not None)

    logger.info('We have split the chembl entries into chunks of ' +
            str(len(ch_het1)) + ', ' + str(len(ch_het2)) + ' and ' +
            str(len(ch_het3)) + ', for easier processing.')

    # filter the struc_map to get rid of domains that are not 
    # interacting with the drug
    chembl_filt_map1 = run_or_pickle("7_chembl_filt_map1", 
                                      drug_targ_res_filter, 
                                      ch_het1, chembl_struct_map)

    # logger.info(len(chembl_filt_map1))

    chembl_filt_map2 = run_or_pickle("7_chembl_filt_map2", 
                                      drug_targ_res_filter, 
                                      ch_het2, chembl_struct_map)
    # logger.info(len(chembl_filt_map2))

    chembl_filt_map3 = run_or_pickle("7_chembl_filt_map3", 
                                      drug_targ_res_filter, 
                                      ch_het3, chembl_struct_map)

    # logger.info(len(chembl_filt_map3))

    chembl_filt_map = dict(chembl_filt_map1.items() + 
                      chembl_filt_map2.items() + chembl_filt_map3.items())

    # logger.info(len(chembl_filt_map))

    # logger.info('The refined dictionary for ChEMBL is: ' +
    #             str(chembl_filt_map))
    
    schis_targ = []
    for drug in chembl_filt_map:
      # logger.info(chembl_filt_map[drug])
      for tar in chembl_filt_map[drug]:
        for ar in chembl_filt_map[drug][tar]:
          for sch in chembl_filt_map[drug][tar][ar]:
            schis_targ.append(sch)
    
    schis_targ = list(set(schis_targ))
    
    logger.info('At this stage we have obtained ' + 
                str(len(chembl_filt_map)) + ' drugs, mapped to '
                + str(len(schis_targ)) +
                ' ' + species_string + ' targets.')


    # percent_identity_chembl = run_or_pickle('7_chembl_percent_identity', 
    #                                       percent_identity, chembl_filt_map)




    # quick fix for drugbank
    drugbank_repo_map2 = AutoVivification()
    for db_id in drugbank_repo_map:
      # logger.info(db_id)
      db_strip = db_id.strip(' ')
      # logger.info(db_strip)
      drugbank_repo_map2[db_strip] = drugbank_repo_map[db_id]



    # for item in drugbank_repo_map:
    #   logger.info(item)
    drugbank_struct_map, drugbank_het_map = (
                  struct_maps(drugbank_repo_map2, drugbank_cluster,
                    uniprot_pdb_w_lig, pdb_cc_dic))
    
    logger.info(len(drugbank_cluster))
    # for dr in drugbank_cluster:
    #   if dr not in drugbank_struct_map:
    #     logger.info(dr)
    logger.info(str(drugbank_cluster['DB07791']))
    logger.info(len(drugbank_struct_map))
    logger.info(len(drugbank_het_map))
    logger.info('check why these two dont match with prev!!!')
    # logger.info(drugbank_struct_map)


    # items1,items2,items3 = \
    #                 zip(*izip_longest(*[iter(drugbank_het_map.items())]*5))
    # db_het1 = dict(item for item in items1 if item is not None)
    # db_het2 = dict(item for item in items2 if item is not None)
    # db_het3 = dict(item for item in items3 if item is not None)
    # db_het4 = dict(item for item in items4 if item is not None)
    # db_het5 = dict(item for item in items5 if item is not None)
    # db_het6 = dict(item for item in items6 if item is not None)
    # db_het7 = dict(item for item in items7 if item is not None)
    # db_het8 = dict(item for item in items8 if item is not None)
    # db_het9 = dict(item for item in items9 if item is not None)
    # db_het10 = dict(item for item in items10 if item is not None)


    # logger.info('We have split the chembl entries into chunks of ' +
    #         str(len(ch_het1)) + ', ' + str(len(ch_het2)) + ' and ' +
    #         str(len(ch_het3)) + ', for easier processing.')

    # # filter the struc_map to get rid of domains that are not 
    # # interacting with the drug
    # db_filt_map1 = run_or_pickle("8_db_filt_map1", 
    #                                   drug_targ_res_filter, 
    #                                   db_het1, drugbank_struct_map)

    # # logger.info(len(chembl_filt_map1))

    # db_filt_map2 = run_or_pickle("8_db_filt_map2", 
    #                                   drug_targ_res_filter, 
    #                                   db_het2, drugbank_struct_map)
    # # logger.info(len(chembl_filt_map2))

    # db_filt_map3 = run_or_pickle("8_db_filt_map3", 
    #                                   drug_targ_res_filter, 
    #                                   db_het3, cdrugbank_struct_map)

    # # logger.info(len(chembl_filt_map3))

    # db_filt_map4 = run_or_pickle("8_db_filt_map4", 
    #                                   drug_targ_res_filter, 
    #                                   db_het4, cdrugbank_struct_map)

    # db_filt_map5 = run_or_pickle("8_db_filt_map5", 
    #                                   drug_targ_res_filter, 
    #                                   db_het5, cdrugbank_struct_map)
    # db_filt_map6 = run_or_pickle("8_db_filt_map6", 
    #                                   drug_targ_res_filter, 
    #                                   db_het6, cdrugbank_struct_map)

    # db_filt_map7 = run_or_pickle("8_db_filt_map7", 
    #                               drug_targ_res_filter, 
    #                               db_het7, cdrugbank_struct_map)

    # db_filt_map8 = run_or_pickle("8_db_filt_map8", 
    #                               drug_targ_res_filter, 
    #                           db_het8, cdrugbank_struct_map)

    # db_filt_map9 = run_or_pickle("8_db_filt_map9", 
    #                           drug_targ_res_filter, 
    #                       db_het9, cdrugbank_struct_map)


    # db_filt_map10 = run_or_pickle("8_db_filt_map10", 
    #                           drug_targ_res_filter, 
    #                       db_het10, cdrugbank_struct_map)

    # chembl_filt_map = dict(chembl_filt_map1.items() + 
    #                   chembl_filt_map2.items() + chembl_filt_map3.items())





    logger.info('----------------------------------------------------------')

  else:
    logger.info('We have skipped STEP ' + str(step))


  ########################################
  ### STEP 9 REPOSITIONING CANDIDATE   ###
  ########################################
  step = step + 1

  if c.steps > (step-1):
    logger.info('------------------------- STEP ' + str(step) + ' ' +
                '-------------------------')

    # DRUG REPOSITIONING CANDIDATE!!
    # info retrieval and alignment
    logger.info('We wish to investigate the repositioning candidate ' + 
                c.repo_candidate + '.')
    # chembl id format
    chembl_format = re.compile('CHEMBL.*')
    # drugbank id format
    drugbank_format = re.compile('DB.*')



    # CHEMBL
    if chembl_format.match(c.repo_candidate):
      # display chembl_struct_map, map with only structural 
      struct_map = chembl_struct_map[c.repo_candidate]
      # logger.info(chembl_het_map[c.repo_candidate])
      het_map = chembl_het_map[c.repo_candidate]

      filt_map = chembl_filt_map[c.repo_candidate]

    # DRUGBANK
    elif drugbank_format.match(c.repo_candidate):
      # display chembl_struct_map, map with only structural 
      struct_map = drugbank_struct_map[c.repo_candidate]
      het_map = drugbank_het_map[c.repo_candidate]

  
    logger.info('The mapping dictionary for the drug is ' + 
                 str(struct_map))

    logger.info('The het dictionary for the drug is ' + 
                 str(het_map))

    logger.info('The filtered dictionary for the drug is ' + 
               str(filt_map))


    # for targ in filt_map:
    #   logger.info(targ)
    #   for arc in filt_map[targ]:
    #     logger.info(arc)

    # fooz = []
    # for drug in chembl_filt_map:
    #   if drug == 'CHEMBL941' or drug == 'CHEMBL1421':
    #     for targ in chembl_filt_map[drug]:
    #       for arch in chembl_filt_map[drug][targ]:
    #         for schis in chembl_filt_map[drug][targ][arc]:
    #           # logger.info(schis)
    #           fooz.append(schis)

    # fooz = list(set(fooz))
    # logger.info(len(fooz))
    # isit = 'G4VP53' in fooz
    # logger.info(isit)
    # # logger.info(len(full_map.values()))


    # # targets we are interested in
    # logger.info('The drug was mapped to ' + str(len(full_map)) + 
    #             ' drug target(s), ' + str(full_map.keys()))

    # logger.info('Of these targets, we are only interested in ' +
    #             str(len(lucky_uniprot)) + 
    #             ' - the one(s) that have crystal structure(s) in complex' + 
    #             ' with ' + str(ref_het) + ', the het group(s) ' +
    #             'the drug is associated to in the clustering.')

    # logger.info('The target(s), associated with their pdb ids, are: ' +
    #             str(lucky_uniprot))


    # # drug target to focus on
    # target_to_align = lucky_uniprot.keys()[c.repo_target_no]

    # # if there is more than oen target, inform on which one we are 
    # focusing on
    # # can change the target number in config file, default is the first one
    # if len(lucky_uniprot) > 1:
    #   logger.info('We are focusing on target number ' + 
    #               str(c.repo_target_no + 1) + ', ' + 
    #               target_to_align + '.')


    # # get list of uniprots, the drug targets and all the ones mapped to it
    # uniprot_to_align = []
    # # add all targets
    # for arch in partial_map[target_to_align]:
    #   for uni in partial_map[target_to_align][arch]:
    #     uniprot_to_align.append(uni)
    #   # logger.info(partial_map[target_to_align][arch])
    
    # #get rid of duplicates
    # uniprot_to_align = list(set(uniprot_to_align))

    # # logger.info(len(uniprot_to_align))
    # logger.info('We are aligning it to all the ' + species_string +
    #             ' targets that were found to be related to it, (' +
    #             str(len(uniprot_to_align)) + ' unique targets).')
    
    # # add the drug target at the beginning of the list
    # uniprot_to_align.insert(0, target_to_align) 
    # # logger.info(uniprot_to_align)

    # #write fasta file from list of uniprots, simply returns name of the file
    # alignment_name = run_or_pickle('9_align', uniprot_to_fasta, 
    #                                 uniprot_to_align)
    
    # logger.info('We have written a fasta file, ' + alignment_name +
    #             ', with the protein sequences not aligned.')


    # # align list of uniprot with t-coffee
    # run_tcoffee(str(alignment_name))

  
    logger.info('----------------------------------------------------------')

  else:
    logger.info('We have skipped STEP ' + str(step))

  ####################################
  ### STEP 10 HOMOLOGY MODELLING   ###
  ####################################
  step = step + 1

  if c.steps > (step-1):
    logger.info('------------------------- STEP ' + str(step) + ' ' +
                '-------------------------')

    logger.info('We wish to build ' + str(c.model_no) +
                ' homology model of the sequence ' + c.model_seq + 
                ', using the template ' + c.model_xray + 
                ' and the alignment file ' + c.model_align + '.')
    


    # run modeller
    # requires alignment file and pdb file in the working dir
    # return sorted list of models
    model_foo = run_or_pickle('10_model_foo', run_modeller, c.model_no,
                              c.model_align, c.model_xray, c.model_seq)

    for i in range(0,len(model_foo)):
      logger.info(str(model_foo[i]['num']) + ' ' +
                  str(model_foo[i]['name']) + str(model_foo[i]['DOPE score']))


    logger.info('We have successfully built ' + 
                str(len(model_foo)) + ' homology model(s).')
    logger.info('The top scoring model is ' + str(model_foo[0]['name']) +
                ', with a DOPE score of ' + str(model_foo[0]['DOPE score']) +
                '.')
  
    logger.info('----------------------------------------------------------')

  else:
    logger.info('We have skipped STEP ' + str(step))
 



  end_time = datetime.now()

  logger.info('The total runtime of the script is: ' + 
              str(end_time - start_time))
  
  logger.info('******************** END OF SCRIPT ' +
              '********************')

############################################################################




############################################################################
### call main function, prevent execution on import
############################################################################
if __name__ == "__main__":
  main()
############################################################################