'''
NAME
    Einfo

VERSION
    1.0

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/Einfo.py

DESCRIPTION
    Provides more information about the data bases

CATEGORY
    Entrez

USAGE
    python Einfo.py

ARGUMENTS
    None

SEE ALSO
    

'''

# Import libraries
from Bio import Entrez
from pprint import pprint

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Use Einfo and read the content
handle = Entrez.einfo()
record = Entrez.read(handle)
handle.close()

# Print
pprint(record["DbList"])
