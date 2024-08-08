'''
NAME
    Esearch_retmax

VERSION
    1.0

TYPE
    Example

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/Esearch_retmax.py

DESCRIPTION
    Provides the number of results in a specific data base and the 
    whole ID list

CATEGORY
    Entrez

USAGE
    python src/Esearch_retmax.py

ARGUMENTS
    None

SEE ALSO
    Esearch
    Esearch_manyfields
    Esearch_author

'''

# Import libraries
from Bio import Entrez

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Search the term in the data base
handle = Entrez.esearch(db = "pubmed", 
                        term = "biopython")
record = Entrez.read(handle)
handle.close()

# Specify to get all the IDs
count = int(record["Count"])
handle = Entrez.esearch(db = "pubmed", 
                        term = "biopython", 
                        retmax = count)
record = Entrez.read(handle)
handle.close()

# Print results
print(len(record["IdList"]))
print(record["IdList"])
