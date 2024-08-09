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
    python src/NCBI_search.py

ARGUMENTS
    

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
                    default="results/esearch_results.txt")

# Store the arguments
args = parser.parse_args()

# Register the database and the terms
data_base = args.database
title = args.title
author = args.author
filename = "results/" + args.output

# Verify that there is enough information to search
try:
    if not title and not author:
        raise InvalidInputError("\nThere is not information to search.")

except InvalidInputError as invalid_input_error:
    print(invalid_input_error.args[0] 
          + " Please, enter at least an author or " 
          + "a term from the title.\n")

else:

    # Verify that the database exists
    try:
        handle = Entrez.einfo()
        record = Entrez.read(handle)
        handle.close()
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

        # Search the terms in the specified database
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

        # Save the IDs in a new file
        with Entrez.esearch(db = data_base, term = term_to_search) as handle:
            record = Entrez.read(handle)

        with open(filename, "w") as file:
            for id in record["IdList"]:
                file.write(f"{id}\n")

        # Let the user know that the search is finished
        print(f"\nThe results are now available in {filename}\n")
