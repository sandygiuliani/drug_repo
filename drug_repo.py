# Copyright 2014 Sandra Giuliani 
# drug_repo.py

# Mapping of known drugs from ChEMBL and DrugBank to their targets/domain 
# architecture to identify suitable drug repositioning candidates for 
# schistosomiasis.
# Please see README.md for more info.




############################################################################
### import modules and define constant variables
############################################################################

# import modules
import sys, re, string, os, fnmatch, shutil

# set up log
import logging

# define CHEMBL_INPUT as the chembldrugs file
CHEMBL_INPUT = 'chembldrugs.txt'
# define CHEMBL_UNIPROT as the chemblID/uniprot mapping file
CHEMBL_UNIPROT = 'chembl_uniprot_mapping.txt'

############################################################################



############################################################################
### HEADER_TAB_COUNT HELPER FUNCTION
############################################################################
# find specific header in first line of tab-separated file 
# and return column number of the header
def header_tab_count(tab_file, header):
  '''(str)->int
  Read first line of tab separated file tab_file and return column number 
  of header
  >>>column_count("DEVELOPMENT_PHASE")
  3
  '''

  # read just first line of tab_file
  with open(tab_file, 'r') as f:
    first = f.readline()
  # set counter to 0
  col_number = 0
  # loop until word matches the header
  while not first.split("\t")[col_number] == header:
    col_number = col_number + 1
  # return column number
  return col_number

############################################################################




############################################################################
### PROCESS_CHEMBL FUNCTION
############################################################################

def processchembl():
  '''read chembl drug input file and return information on number of drugs 
  and headers'''
  # open chembldrugs.txt for reading
  input_file = open(CHEMBL_INPUT, 'r') 
  lines = input_file.readlines()
  input_file.close()
  logging.info('The number of drugs listed in the input file is '
              + str(len(lines)-1))
  # get the headers
  headers = lines[0]
  # remove duplicate spaces
  headersnospace = " ".join(headers.split())
  # print headers list in lowercase
  # print('The headers are: '+headersnospace.lower())

  col_phase = header_tab_count(CHEMBL_INPUT,"DEVELOPMENT_PHASE")
  col_type = header_tab_count(CHEMBL_INPUT,"DRUG_TYPE")
  
  logging.info('The column with the development_phase info is the '+ str(col_phase+1)
       + 'th; the column with the drug_type info is the '
       + str(col_type+1) + 'th.')

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

  logging.info('Number of drugs in phase 1 is: '+ str(phase1)+'; in phase 2 is: ' 
        +str(phase2)+'; in phase 3 is: ' +str(phase3)+'; in phase 4 is: ' 
        +str(phase4)+'; in unknown phase is: ' +str(phase_unknown)+'.')

  # open the file to write to
  stripped = open('chembldrugs_stripped.txt', 'w')
  # set counter to zero, this is just to know how many lines we end up writing
  total_stripped_lines = 0
  # look over lines, excluding header
  for y in range(1,len(lines)):
    # tab separate each row
    rowsplit2 = lines[y].split("\t")
    # check if they are phase 4 or unknown
    # can also add phase 3 from here rowsplit2[col_phase] == '3')
    if (rowsplit2[col_phase] == '4') or (rowsplit2[col_phase] == ''): 
      # increase the useless counter
      total_stripped_lines = total_stripped_lines + 1
      # write to the file the stripped lines
      stripped.write(lines[y])
      # print(rowsplit2[column])
  # print friendly statement
  logging.info('We have written to the file chembl_stripped.txt' + 
      'only the entries in phase 3, 4 or with unkown phase' +
      ', for a total number of '+ str(total_stripped_lines)+' drugs.')
  # close the file we wrote to
  stripped.close()
  
  # we open the stripped file for reading
  stripped2 = open('chembldrugs_stripped.txt', 'r')
  # reading lines
  lines2 = stripped2.readlines()
  # closing the stripped file
  stripped2.close()
  typecount = 0
  # look over, note here there is no header
  for x in range(len(lines2)):
    # tab separate
    rowsplit3 = lines2[x].split("\t")
    # drug_type = ''
    if not rowsplit3[col_type] == 'Synthetic Small Molecule':
      typecount = typecount +1
      # print(rowsplit3[col_type])
  logging.debug(typecount)
  
############################################################################




############################################################################
### MAIN FUNCTION
############################################################################

# call processchembl, 
# TO ADD: call processdrugbank, merge uniprot codes

def main():
  logging.basicConfig(filename='log_drug_repo.log', level=logging.DEBUG)
  processchembl()

############################################################################




############################################################################
### call main function, prevent excecution on import
############################################################################
if __name__ == "__main__":
  main()
############################################################################
