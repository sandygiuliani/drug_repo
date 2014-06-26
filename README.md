# drug_repo #
_a repository of Python scripts for drug repositioning_


Find here Python scripts for a drug repositiong bio-/chemoinformatics project
applied to schistosomiasis.  



**FAQ**  
Q. What is drug repositioning?   
A. The usage of a known drug for a different therapeutic indication. If you are not familiar with this at all, try [Wikipedia](http://en.wikipedia.org/wiki/Drug_repositioning)  
Q. What is schistosomiasis?  
A. A very nasty parasitic disease affecting over 200 million people. Learn more about schistosomiasis on the [World Health Organization website](http://www.who.int/topics/schistosomiasis/en/)  
Q. How does the tool work?  
A. By mapping! known drugs -> their targets -> their domain architecture -> parasite targets  



**CONTENTS**  
* drug_repo.py - Python script that reads input files (chemb/drugbank), filters data, extracts relevant info for mapping with domain architecture info. It is being developed at the moment.
* chembl\_drugs.txt - ChEMBL drugs; downloaded from [ChEMBL](http://www.ebi.ac.uk/chembl/drugstore/)
* chembl\_drugtargets.txt - ChEMBL drug targets; downloaded from [ChEMBL](http://www.ebi.ac.uk/chembl/drug/targets/) and manually edited to strip a newline character at lines 383/384.  
* chembl\_uniprot\_mapping.txt - ChEMBL uniprot mapping, chembl ID to UniProt codes; downloaded from the ChEMBL 18 release page: ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_18/  
* all\_target\_ids\_all.csv - DrugBank Drug Target Identifiers/All Drugs; downloaded from [DrugBank](http://www.drugbank.ca/downloads#protein-identifiers)  
* small\_molecule\_target\_ids\_all.csv - DrugBank Drug Target Identifier/Small Molecule Drugs; downloaded from [DrugBank](http://www.drugbank.ca/downloads#protein-identifiers)  
* uniprot_pdb.*sv (csv and tsv) - Uniprot to pdb mapping file; downloaded from [SIFTS](http://www.ebi.ac.uk/pdbe/docs/sifts/quick.html), accessed 12/06/2014  
* het_pairs.lst - pdb to Het groups mapping file; downloaded from [PDBsum downloads](http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?doc=TRUE&template=downloads.html&pdbcode=n/a), accessed 12/06/2014  
* lig_pairs.lst - pdb to ligand mapping file; downloaded from [PDBsum downloads](http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?doc=TRUE&template=downloads.html&pdbcode=n/a), accessed 17/06/2014  
* Components-smiles-*.smi - chemical components dictionary in smiles format (stereo for those including stereocenters, oa for OpenEye and cactvs for CACTVS-generated); downloaded from [RCSB Ligand Expo Downloads](http://ligand-expo.rcsb.org/ld-download.html), in the SMILES/InChi data files, accessed 20/06/2014  
* \*.p - any .p file is a pickle (for caching purposes) if there are any around, it means I need to move them across machines, please ignore  



**REQUIREMENTS**  
* ArchIndex/ArchSchema - kindly provided by Dr Laskowski. For more information, please visit the [ArchSchema website](http://www.ebi.ac.uk/thornton-srv/databases/archschema), or read the [main reference for ArchSchema](http://www.ncbi.nlm.nih.gov/pubmed/20299327)  
* SMSD (Small Molecule Subgraph Detector). For more information, please visit the [SMSD website](http://www.ebi.ac.uk/thornton-srv/software/SMSD/), the [GitHub repository](https://github.com/asad/SMSD), or read the [main reference for SMSD](http://www.jcheminf.com/content/1/1/12)  



**CHECKLIST**  
drug_repo.py development:
- [x] ChEMBL processing - clinical phase filter 
- [x] ChEMBL processing - small molecule filter 
- [x] ChEMBL processing - uniprot mapping
- [x] DrugBank data processing
- [x] write everything to log
- [x] ChEMBL/DrugBank data merge
- [x] domain architecture mapping (ArchSchema)
- [x] schistosoma domains
- [x] ChEMBL map
- [x] DrugBank map
- [ ] structure (pdb) filter
- [ ] targets shortlist



**LICENSE**  
Copyright &copy; 2014 Sandra Giuliani  
This repository is licensed under the terms of the MIT license. Please see LICENSE.md for more information.  
The MIT license is approved by the [Open Source Initiative](http://opensource.org/licenses)  



**DISCLAIMER**  
THIS IS A WORK IN PROGRESS. I AM NEW TO GITHUB AND NEW TO PROGRAMMING IN GENERAL! Feedback is welcome but please be kind.  


**CONTACT**  
Drop me a line at: sandraxgiuliani [at] gmail [dot] com  
You might also want to follow me on Twitter [@sandygiuliani](https://twitter.com/sandygiuliani) or visit my [personal website](http://www.sandragiuliani.com/).  
