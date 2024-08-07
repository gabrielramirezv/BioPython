'''
NAME
    amino_acid_content

VERSION
    1.0

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/Programacion_en_Python/blob/master/src/amino_acid_content.py

DESCRIPTION
    Provides more information about the data bases

CATEGORY
    Entrez

USAGE
    python Einfo.py

ARGUMENTS
    None

SEE ALSO
    AT_GC_percentage

'''

# Import libraries
from Bio import Entrez
from pprint import pprint

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Use Einfo and read the content
handle = Entrez.einfo()
result = handle.read()
handle.close()

# Print
pprint(result)
