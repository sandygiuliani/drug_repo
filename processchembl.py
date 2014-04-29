#Copyright 2014 Sandra Giuliani 
#processchembl reads and extracts
#information from chembl drug file chembldrugs.txt make sure the input file
#matches the one downloaded from chembl, from the drug tab the txt file is tab
#separated and has a first line of headers the header with the drug
#development phase information (phase 1,2,3 or 4) must be called
#'DEVELOPMENT_PHASE' the column with the development_phase info should be the
#4th (3rd index), but in case they change it in future releases, the program
#makes no assumption about the number and searches for the string.

#tested with chembl_18 release (April 2014) 
#the total number of drugs listed is 10,406


#import modules
import sys, re, string, os, fnmatch, shutil

#define constant variable CHEMBL_INPUT for chembldrugs input
CHEMBL_INPUT = 'chembldrugs.txt'

###
#PROCESS CHEMBL FUNCTION CALL
###

def processchembl():
  '''reads chembl drug input file and returns information on number of drugs 
  and headers'''
  #opens chembldrugs.txt for reading
  input_file = open(CHEMBL_INPUT, 'r') 
  lines = input_file.readlines()
  input_file.close()
  print('The number of drugs listed in the input file is '+str(len(lines)-1))
  #gets the headers
  headers = lines[0]
  #removes duplicate spaces
  headersnospace = " ".join(headers.split())
  #prints headers list in lowercase
  #print('The headers are: '+headersnospace.lower())


  #CLINICAL PHASE FILTER
  
  #splits header row according to tab separators
  headersplit = lines[0].split("\t")
  #print(headersplit[3] == "DEVELOPMENT_PHASE")
  col_phase = 0
  while not (headersplit[col_phase] == "DEVELOPMENT_PHASE"):
    col_phase = col_phase +1
  col_type = 0
  while not (headersplit[col_type] == "DRUG_TYPE"):
    col_type = col_type +1
  #the value col_phase now refers to the index of the development_phase,
  #the value col_type to the drug_type
  #print out which column we are referring to 
  print('The column with the development_phase info is the '+str(col_phase+1)
      +'th; the column with the drug_type info is the '+str(col_type+1)+'th.')
  

  #CLINICAL PHASE FILTER

  ##this bit counts the drugs in the 4+1 classes of clinical phase
  #set counters for clinical phases to zero
  phase1 = 0
  phase2 = 0
  phase3 = 0
  phase4 = 0
  phase_unknown =0
  #iterate over rows, excluding the header row
  for i in range(1,len(lines)):
    #tab separate each row
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

  print('Number of drugs in phase 1 is: '+ str(phase1)+'; in phase 2 is: ' 
        +str(phase2)+'; in phase 3 is: ' +str(phase3)+'; in phase 4 is: ' 
        +str(phase4)+'; in unknown phase is: ' +str(phase_unknown)+'.')

  #open the file to write to
  stripped = open('chembldrugs_stripped.txt', 'w')
  #set counter to zero, this is just to know how many lines we end up writing
  total_stripped_lines = 0
  #look over lines, excluding header
  for y in range(1,len(lines)):
    #tab separate each row
    rowsplit2 = lines[y].split("\t")
    #check if they are phase 4 or unknown
    #can also add phase 3 from here rowsplit2[col_phase] == '3')
    if (rowsplit2[col_phase] == '4') or (rowsplit2[col_phase] == ''): 
      #increase the useless counter
      total_stripped_lines = total_stripped_lines + 1
      #write to the file the stripped lines
      stripped.write(lines[y])
      #print(rowsplit2[column])
  #print friendly statement
  print('We have written to the file chembl_stripped.txt' + 
      'only the entries in phase 3, 4 or with unkown phase' +
      ', for a total number of '+ str(total_stripped_lines)+' drugs.')
  #close the file we wrote to
  stripped.close()
  
  #we open the stripped file for reading
  stripped2 = open('chembldrugs_stripped.txt', 'r')
  #reading lines
  lines2 = stripped2.readlines()
  #closing the stripped file
  stripped2.close()
  typecount = 0
  #look over, note here there is no header
  for x in range(len(lines2)):
    #tab separate
    rowsplit3 = lines2[x].split("\t")
    #drug_type = ''
    if not rowsplit3[col_type] == 'Synthetic Small Molecule':
      typecount = typecount +1
      #print(rowsplit3[col_type])
  print(typecount)
  

###
#END OF MAIN FUNCTION
###




#calls processchembl, prevents excecution on import
if __name__ == "__main__":
  processchembl()
