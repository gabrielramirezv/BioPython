'''
NAME
    Uniprot_accession

VERSION
    2.0

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/Uniprot_accession.py

DESCRIPTION
    Searches a specific protein in different databases, and it returns 
    the IDs and accessions.

CATEGORY
    ExPASy
    
REQUIREMENTS
    - Software - 
        You must have Python and BioPython installed
    - Repository structure -
        For using this script successfully, you should have the 
        following repository structure:
        |- data      # It contains data to be used as input
        |- docs      # It contains documentation
        |- results   # Output files will be created in this directory
        \- src       # The script must be executed from this directory

USAGE
    python src/Uniprot_accession.py

ARGUMENTS
                        
VARIABLES DICTIONARY
    

SEE ALSO
    Esearch

'''

# Import libraries
import os
from Bio import Entrez, SeqIO, ExPASy

# Register an email
Entrez.email = "gramirez@lcg.unam.mx"

# Search in the protein database
with Entrez.esearch(db="protein", 
                   term="DEFA[gene] AND Aedes aegypti[Orgn]") as handle:
    # Read the results
    record = Entrez.read(handle)

# Get ID
prot_id = record["IdList"]

# Extract a GenBank from the protein database
with Entrez.efetch(db="protein", 
                   id=prot_id, 
                   rettype="gb", 
                   retmode="text") as handle:
    # Read the results
    record = SeqIO.parse(handle, "genbank")
    for protein in record:
        if "UniProtKB" in protein.annotations['db_source']:
            defa_prot = protein.annotations['accessions']

# Create the files with the results and inform the user
if defa_prot:
    print("\nThe following files have been created in the directory results/\n")
    for accession in defa_prot:
        with open(f"results/{accession}.txt", 'w') as accession_file:
            handle = ExPASy.get_sprot_raw(accession)
            accession_file.write(handle.read())
        print(f"{accession}.txt")
else:
    print("\nNo results found.\n")

