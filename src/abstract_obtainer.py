'''
NAME
    abstract_obtainer

VERSION
    2.0

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/abstract_obtainer.py

DESCRIPTION
    Receives a file with IDs of articles of interest and creates a file
    with the abstracts of those articles, including IDs of other 
    papers that cited them.

CATEGORY
    Entrez
    
REQUIREMENTS
    - Software - 
        You must have Python and BioPython installed
    - Repository structure -
        For using this script successfully, you should have the 
        following repository structure:
        |- data      # It contains data to be used as input
        |- docs      # It contains documentation
        |- results   # Output files will be created in this directory
        \- src       # The script must be executed from this directory

USAGE
    python src/abstract_obtainer.py [-h] [-e EMAIL] [-i INPUTFILE] 
        [-o OUTPUTFILE]

ARGUMENTS
    -h, --help          Show a help message and exit
    -e EMAIL, --email EMAIL
                        Email to use in Entrez
    -i INPUTFILE, --inputfile INPUTFILE
                        Path to the input file
    -o OUTPUTFILE, --outputfile OUTPUTFILE
                        Path to the output file
                        
VARIABLES DICTIONARY
    

SEE ALSO
    Esearch
    Elink

'''

# Import libraries
import argparse
from Bio import Entrez

# Create parser
parser = argparse.ArgumentParser(description="Receives a file with " 
                                 + "IDs of articles of interest and " 
                                 + "creates a file with the abstracts "
                                 + "of those articles, including IDs "
                                 + "of other papers that cited them.")

# Store the arguments
parser.add_argument("-e", "--email",
                    help="Email to use in Entrez",
                    type=str,
                    required=False,
                    default="gramirez@lcg.unam.mx")

parser.add_argument("-i", "--inputfile",
                    help="Path to the input file",
                    type=str,
                    required=False,
                    default="results/j_collado_papers_about_regulation.txt")

parser.add_argument("-o", "--outputfile",
                    help="Path to the output file",
                    type=str,
                    required=False,
                    default="results/abstracts.txt")

args = parser.parse_args()
email = args.email
input_file = args.inputfile
output_file = args.outputfile

# Register an email account
Entrez.email = email

# Take the file with the article IDs as the input file
with open(input_file, "r") as ids_list:
    articles_ids = ids_list.read()
    articles_ids = articles_ids.replace(r"\n", ",")
try:
    articles_ids = articles_ids.split()
    # Save the abstract of the articles
    with open(output_file, 'w') as out_handle:
        article_number = 0
        for article in articles_ids:
            article_number += 1
            with Entrez.efetch(db="pubmed",
                            id=article,
                            rettype = "abstract",
                            retmode="text") as fetch_handle:
                abstracts = fetch_handle.read().replace("1. ", 
                                                        f"{article_number}. ")
                
                # Save the IDs of other articles that cite it
                results = Entrez.read(Entrez.elink(dbfrom="pubmed", 
                                                   db="all",
                                                   LinkName="pubmed_pmc_refs", 
                                                   from_uid=article))
                
                # Write the abstracts with the link IDs
                try:

                    pmc_ids = [link["Id"] for link 
                               in results[0]["LinkSetDb"][0]["Link"]]
                    out_handle.write(f"{abstracts}\nCites: " 
                                     + f"{', '.join(pmc_ids)}\n\n\n\n")
                except:
                    out_handle.write(f"{abstracts}\nCites: There are no " + 
                                     "cites to this article.\n\n\n\n")

    # Let the user know that the output file has been created
    print("\nThe file results/abstracts.txt has been created successfully.\n")
except:
    # If there was a problem writing the file, tell the user what 
    # article could not be processed
    print("\n An error ocurred while creating the file. There was a " + 
          f"problem processing article number {article_number}.\n")
    