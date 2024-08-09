'''
NAME
    Efetch

VERSION
    1.0

TYPE
    Example

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/Efetch.py

DESCRIPTION
    Gets records in the specified format from the database

CATEGORY
    Entrez

USAGE
    python src/Efetch.py

ARGUMENTS
    None

SEE ALSO
    Einfo
    Esearch
    Esummary
    Espell

'''

# Import libraries
from Bio import Entrez, SeqIO

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Search the sequence in the database
handle = Entrez.efetch(db = "nucleotide", 
                       id = "HE805982",
                       rettype = "gb",
                       retmode = "text")
record = SeqIO.read(handle, "genbank")
handle.close()

# Print results
print(record)
