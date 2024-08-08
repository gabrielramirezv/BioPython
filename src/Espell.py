'''
NAME
    Espell

VERSION
    1.0

TYPE
    Example

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/Espell.py

DESCRIPTION
    Corrects the spelling of a term to search

CATEGORY
    Entrez

USAGE
    python src/Espell.py

ARGUMENTS
    None

SEE ALSO
    Einfo
    Esearch
    Esummary
    Efetch

'''

# Import libraries
from Bio import Entrez

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Enter the term to search
handle = Entrez.espell(term="biopythooon")
record = Entrez.read(handle)
handle.close()

# Imprimir resultados
print(record["Query"])
print(record["CorrectedQuery"])
