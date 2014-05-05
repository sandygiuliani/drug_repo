# drug_repo #
_a repository of Python scripts for drug repositioning_


Find here Python scripts for a drug repositiong bio-/chemoinformatics project
applied to schistosomiasis.  
What is drug repositioning? The usage of a known drug for a different therapeutic indication. For a few examples, see http://en.wikipedia.org/wiki/Drug_repositioning  
What is schistosomiasis? A very nasty parasitic disease affecting over 200 million people. Learn more at: http://www.who.int/topics/schistosomiasis/en/

**CONTENTS**
* drug_repo.py - Python script that reads input files (chemb/drugbank), filters data, extracts relevant info for mapping with domain architecture info. It is being developed at the moment.
* chembl\_drugs.txt - ChEMBL drugs. From www.ebi.ac.uk/chembl/drugstore ('browse drugs')
* chembl\_drugtargets.txt - ChEMBL drug targets. From www.ebi.ac.uk/chembl/drug/targets ('browse drug targets')  
* chembl\_uniprot\_mapping.txt - ChEMBL uniprot mapping, chembl ID to UniProt codes. From ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_18/

**CHECKLIST**  
drug_repo.py development:
- [x] ChEMBL processing - clinical phase filter 
- [x] ChEMBL processing - small molecule filter 
- [x] ChEMBL processing - uniprot mapping
- [ ] DrugBank data processing
- [x] write everything to log instead of .txt output
- [ ] ChEMBL/DrugBank data merge
- [x] domain architecture mapping (ArchSchema)
- [x] schistosoma domains

**DISCLAIMER:**  
THIS IS A WORK IN PROGRESS. I AM NEW TO GITHUB AND NEW TO PROGRAMMING IN GENERAL! Feedback welcome but please be kind.

Contact me at: sandraxgiuliani@gmail.com

**LICENSE**  
Copyright &copy; 2014 Sandra Giuliani  
This repository is licensed under the terms of the MIT license. Please see LICENSE.md for more information.  
The MIT license is approved by the Open Source Initiative http://opensource.org/licenses


