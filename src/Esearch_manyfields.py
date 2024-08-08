'''
NAME
    Esearch_manyfields

VERSION
    1.0

TYPE
    Example

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/Esearch_manyfields.py

DESCRIPTION
    Provides the count and the list of IDs related with a specific
    term

CATEGORY
    Entrez

USAGE
    python src/Esearch_manyfields.py

ARGUMENTS
    None

SEE ALSO
    Esearch
    Esearch_author
    Esearch_retmax

'''

# Import libraries
from Bio import Entrez

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"
term_to_search = "(Aedes[Title] OR Aedes[All Fields])"
term_to_search += "AND (RNA-Seq[Title] OR transcriptomic[Title] OR "
term_to_search += "transcriptome[Title] OR sequencing[Title])"

# Search in the data base
handle = Entrez.esearch(db = "pubmed", term = term_to_search)
record = Entrez.read(handle)
handle.close()

# Print results
print(record["Count"])
print(record["IdList"])
