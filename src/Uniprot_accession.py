'''
NAME
    Uniprot_accession

VERSION
    3.0

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
    

SEE ALSO
    Esearch

'''

# Import libraries
import os, argparse
from Bio import Entrez, SeqIO, ExPASy

# Create parser
parser = argparse.ArgumentParser(description="Search taxonomy files " 
                                 + "for two different species and " 
                                 + "compare their lineages")

# Store the arguments
parser.add_argument("-e", "--email",
                    help="Email to use in Entrez",
                    type=str,
                    required=False,
                    default="gramirez@lcg.unam.mx")

parser.add_argument("-t", "--term",
                    help="Term to search",
                    type=str,
                    required=False,
                    default="DEFA[gene] AND Aedes aegypti[Orgn]")

args = parser.parse_args()
email = args.email
term_to_search = args.term

# Register an email
Entrez.email = email

# Search in the protein database
with Entrez.esearch(db="protein", 
                   term=term_to_search) as handle:
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
try:
    if defa_prot:
        print("\nThe following files have been created in the " 
              + "directory results/\n")
        for accession in defa_prot:
            with open(f"results/{accession}.txt",
                       'w') as accession_file:
                handle = ExPASy.get_sprot_raw(accession)
                accession_file.write(handle.read())
            print(f"{accession}.txt")
    else:
        print("\nNo results found.\n")
        
# If there is not a results directory, create it
except OSError:
    os.mkdir("results")
    if defa_prot:
        for accession in defa_prot:
            with open(f"results/{accession}.txt",
                       'w') as accession_file:
                handle = ExPASy.get_sprot_raw(accession)
                accession_file.write(handle.read())
            print(f"{accession}.txt")
