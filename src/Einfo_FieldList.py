'''
NAME
    Einfo

VERSION
    2.0

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
    -h, --help          Show a help message and exit
    -db DATABASE, --database DATABASE
                        DataBase that you want to access to

SEE ALSO
    Einfo

'''

# Import libraries
import argparse
from Bio import Entrez

# Define errors
class LackOfInputError(Exception):
    pass

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Create parser
parser = argparse.ArgumentParser(description="Script that access to the "
                                 + "FieldList of a specific data base")

# Store the arguments
parser.add_argument("-db", "--database",
                    help="DataBase that you want to access to",
                    type=str,
                    required=True)

args = parser.parse_args()

# Register the data base
data_base = args.database

# Use Einfo and read the content
handle = Entrez.einfo(db = data_base)
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
