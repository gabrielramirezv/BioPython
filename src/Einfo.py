'''
NAME
    Einfo

VERSION
    3.1

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
    Einfo_FieldList

'''

# Import libraries
from Bio import Entrez

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Use Einfo and read the content
handle = Entrez.einfo()
record = Entrez.read(handle)
handle.close()

# Print
for database in record["DbList"]:
    print(database)
