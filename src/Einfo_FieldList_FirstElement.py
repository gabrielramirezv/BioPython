'''
NAME
    Einfo_FieldList_FirstElement

VERSION
    1.0

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/Einfo_FieldList_FirstElement.py

DESCRIPTION
    Provides more information about the first element of FieldList 
    in DbInfo

CATEGORY
    Entrez

USAGE
    python src/Einfo_FieldList_FirstElement.py [-h] [-db DATABASE]

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
class InvalidInputError(Exception):
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
                    required=False,
                    default="pubmed")

args = parser.parse_args()

# Register the data base
data_base = args.database

# Use Einfo and read the content
handle = Entrez.einfo()
record = Entrez.read(handle)
handle.close()

# Verify that it is a valid data base
try:
    database_is_valid = False
    for database in record["DbList"]:
            if data_base == database:
                database_is_valid = True
                break
    if database_is_valid == False:
        raise InvalidInputError("\nDataBase is not valid.")
except InvalidInputError as invalid_input_error:
    print(invalid_input_error.args[0] 
          + " Please, select one of these data bases:\n")
    for database in record["DbList"]:
        print(database)

else:

    # Use Einfo and read the content of the data base
    handle = Entrez.einfo(db = data_base)
    record = Entrez.read(handle)
    handle.close()

    # Access to FieldList
    fields = record["DbInfo"]["FieldList"]

    # Print the content of the first element
    print("NAME\tFULL NAME\tDESCRIPTION")
    name = fields[0]["Name"]
    full_name = fields[0]["FullName"]
    description = fields[0]["Description"]
    print(f"{name}\t{full_name}\t{description}")
        