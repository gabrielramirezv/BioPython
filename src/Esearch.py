'''
NAME
    Esearch

VERSION
    1.0

TYPE
    Example

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/Esearch.py

DESCRIPTION
    Provides the number of results of your search in the specified 
    data base

CATEGORY
    Entrez

USAGE
    python src/Esearch.py

ARGUMENTS
    None

SEE ALSO
    Einfo

'''

# Import libraries
from Bio import Entrez

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Search the specified term in the data base
handle = Entrez.esearch(db = "pubmed", term = "biopython")
record = Entrez.read(handle)
handle.close()

# Print the results
print(record["Count"])
