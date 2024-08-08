'''
NAME
    Esummary

VERSION
    1.0

TYPE
    Example

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/Esummary.py

DESCRIPTION
    Gets a summary of the information for a list of IDs

CATEGORY
    Entrez

USAGE
    python src/Esummary.py

ARGUMENTS
    None

SEE ALSO
    Einfo
    Esearch
    Espell
    Efetch

'''

# Import libraries
from Bio import Entrez

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Search the ID list in the database
handle = Entrez.esummary(db="taxonomy", id="9913,30521")
record = Entrez.read(handle)
handle.close()

# Print results
print(len(record))
print(record[0]["Id"])