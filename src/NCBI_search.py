'''
NAME
    NCBI_search

VERSION
    1.0

TYPE
    Exercise

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/NCBI_search.py

DESCRIPTION
    Search in a database from NCBI according to the author and terms 
    in the title

CATEGORY
    Entrez

USAGE
    python src/NCBI_search.py [-h] [-db DATABASE] 
    [-a AUTHOR] [-t TITLE WORDS SEPARATED BY COMMAS] [-o OUTPUT]

ARGUMENTS
    -h, --help              show this help message and exit
    -db DATABASE, --database DATABASE
                            DataBase that you want to access to
    -a AUTHOR, --author AUTHOR
                            Term to search in the field "AUTHOR"
    -t TITLE, --title TITLE
                            Terms to search in the title separated by commas
    -o OUTPUT, --output OUTPUT
                            Name of the output file

SEE ALSO
    Esearch

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
parser = argparse.ArgumentParser(description="Search in a database " 
                                 + "from NCBI according to the author" 
                                 + " and terms in the title")

# Define the arguments
parser.add_argument("-db", "--database",
                    help="DataBase that you want to access to",
                    type=str,
                    required=False,
                    default="pubmed")

parser.add_argument("-a", "--author",
                    help="Term to search in the field \"AUTHOR\"",
                    type=str,
                    required=False)

parser.add_argument("-t", "--title",
                    help="Terms to search in the title separated by commas",
                    type=str,
                    required=False)

parser.add_argument("-o", "--output",
                    help="Name of the output file",
                    type=str,
                    required=False,
                    default="esearch_results.txt")

# Store the arguments
args = parser.parse_args()

# Register the database and the terms
data_base = args.database
title = args.title
author = args.author
filename = "results/" + args.output

# Verify that there is at least a word to search
try:
    if not title and not author:
        raise InvalidInputError("\nThere is not information to search.")

# If there are not words to search, ask the user for some
except InvalidInputError as invalid_input_error:
    print(invalid_input_error.args[0] 
          + " Please, enter at least an author or " 
          + "a term from the title.\n")

else:

    # Verify that the database exists
    try:
        # Get information about the available databases
        handle = Entrez.einfo()
        record = Entrez.read(handle)
        handle.close()
        database_is_valid = False
        # Search the database in the existing databases list
        for database in record["DbList"]:
                if data_base == database:
                    database_is_valid = True
                    break
        # If the database does not exist, give to the user a list of
        # the existing databases
        if database_is_valid == False:
            raise InvalidInputError("\nDataBase is not valid.")
        
    except InvalidInputError as invalid_input_error:
        print(invalid_input_error.args[0] 
            + " Please, select one of these data bases:\n")
        for database in record["DbList"]:
            print(database)

    else:

        # Create the term to search, using the title terms, the author,
        # or both if provided by the user
        if title and author:
            terms_in_title = title.replace(",", "[Title] OR ")
            terms_in_title += "[Title]"
            term_to_search = f"{author}[Author] AND ({terms_in_title})"
        elif title and not author:
            terms_in_title = title.replace(",", "[Title] OR ")
            terms_in_title += "[Title]"
            term_to_search = f"{terms_in_title}"
        else:
            term_to_search = author

        # Search the term in the database and get the information
        with Entrez.esearch(db = data_base, 
                            term = term_to_search) as handle:
            record = Entrez.read(handle)

        # Create the file and write the IDs in it
        with open(filename, "w") as file:
            for id in record["IdList"]:
                file.write(f"{id}\n")

        # Let the user know that the search is finished and where
        # the results are
        print(f"\nThe results are now available in {filename}\n")
