'''
NAME
    Esearch_author

VERSION
    1.0

TYPE
    Example

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/Esearch_author.py

DESCRIPTION
    Provides the number of results of a search by author in a specific 
    data bases

CATEGORY
    Entrez

USAGE
    python src/Esearch_author.py

ARGUMENTS
    None

SEE ALSO
    Esearch
    Esearch_manyfields
    Esearch_retmax

'''

# Import libaries
from Bio import Entrez

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Search with stablished parameters
handle = Entrez.esearch(db = "pubmed", 
                        term = "Mateo-Estrada V",
                        field = "AUTH")
record = Entrez.read(handle)
handle.close()

# Print the count of results and the ID list
print(record["Count"])
print(record["IdList"])
