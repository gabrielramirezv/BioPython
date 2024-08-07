'''
NAME
    Einfo

VERSION
    1.0

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/Einfo_FieldList.py

DESCRIPTION
    Provides more information about FieldList in DbInfo

CATEGORY
    Entrez

USAGE
    python Einfo_FieldList.py

ARGUMENTS
    None

SEE ALSO
    Einfo

'''

# Import libraries
from Bio import Entrez

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Use Einfo and read the content
handle = Entrez.einfo(db = "pubmed")
record = Entrez.read(handle)
handle.close()

# Access to FieldList
fields = record["DbInfo"]["FieldList"]

# Print the content
print("NAME\tFULL NAME\tDESCRIPTION")
for field in fields:
    name = field["Name"]
    full_name = field["FullName"]
    description = field["Description"]
    print(f"{name}\t{full_name}\t{description}")
