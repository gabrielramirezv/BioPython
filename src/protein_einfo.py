'''
NAME
    protein_einfo

VERSION
    0.1

TYPE
    Exercise

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/protein_einfo.py

DESCRIPTION
    Gets the description of the field "ECNO" and the link 
    "protein_protein_small_genome" from the protein database

CATEGORY
    Entrez

USAGE
    python src/protein_einfo.py

ARGUMENTS
    None

SEE ALSO
    Einfo

'''

# Import libraries
from Bio import Entrez

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Get the information from the database
handle = Entrez.einfo(db = "protein")
record = Entrez.read(handle)
handle.close()

# Get the description of ECNO
fields = record["DbInfo"]["FieldList"]
for field in fields:
    if field["Name"] == "ECNO":
        ecno_description = field["Description"]
        break

# Get the description of protein_protein_small_genome
fields = record["DbInfo"]["LinkList"]
for field in fields:
    if field["Name"] == "protein_protein_small_genome":
        protein_protein_small_genome_description = field["Description"]
        break

# Print results
print(f"ECNO\t{ecno_description}\n" +
      "protein_protein_small_genome\t" + 
      f"{protein_protein_small_genome_description}")
