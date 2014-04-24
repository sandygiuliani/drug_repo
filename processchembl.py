#Copyright 2014 Sandra Giuliani
#processchembl.py reads and extracts information from chembl drug file chembldruginput.txt
#make sure the input file matches the one downloaded from chembl, from the drug tab
#the txt file is tab separated and has a first line of headers
#in April 2014 the total number of drugs listed is 10,406


#import modules
import sys, re, string, os, fnmatch, shutil

#main function
def processchembl():
  '''(str)->str
  reads chembl drug input file and returns information on number of drugs and headers'''
  #opens chembldruginput.txt for reading
  input_file = open('chembldruginput.txt', 'r') 
  lines = input_file.readlines()
  print('The number of drugs listed in the input file is '+str(len(lines)-1))
  #gets the headers
  headers = lines[0]
  headersnospace = " ".join(headers.split())
  print('The headers are: '+headersnospace.lower())
  #print(headers)

  #here I want to add the drug type and clinical phase filter
  for i in range(len(lines)+1):
    return lines[i]

  input_file.close()



#calls processchembl, prevents excecution on import
if __name__ == "__main__":
  processchembl()
