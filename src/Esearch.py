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
    Provides more information about the data bases

CATEGORY
    Entrez

USAGE
    python Esearch.py

ARGUMENTS
    None

SEE ALSO
    Einfo

'''

# Import libraries
from Bio import Entrez

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Realizar la busqueda
handle = Entrez.esearch(db = "pubmed", term = "biopython")
record = Entrez.read(handle)
handle.close()

# Imprimir el conteo de resultados
print(record["Count"])
