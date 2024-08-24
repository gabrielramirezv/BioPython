'''
NAME
    abstract_obtainer

VERSION
    1.0

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
    python src/abstract_obtainer.py

ARGUMENTS
                        
VARIABLES DICTIONARY
    

SEE ALSO
    Esearch
    Elink

'''

# Import libraries
from Bio import Entrez

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Take the file with the article IDs as the input file
with open("results/j_collado_papers_about_regulation.txt", "r") as ids_list:
    articles_ids = ids_list.read()
    articles_ids = articles_ids.replace(r"\n", ",")
try:
    articles_ids = articles_ids.split()
    # Save the abstract of the articles
    with open("results/abstracts.txt", 'w') as out_handle:
        for article in articles_ids:
            with Entrez.efetch(db="pubmed",
                            id=article,
                            rettype = "abstract",
                            retmode="text") as fetch_handle:
                abstracts = fetch_handle.read()
                # Save the IDs of other articles that cite it
                results= Entrez.read(Entrez.elink(dbfrom="pubmed", db="all",
                                        LinkName="pubmed_pmc_refs", from_uid=article))
                try:
                    pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
                    out_handle.write(f"{abstracts}\nCites: {', '.join(pmc_ids)}\n\n\n\n")
                except:
                    out_handle.write(f"{abstracts}\nCites: This article has not any cites in the available databases.\n\n\n\n")

# Let the user know that the output file has been created
    print("\n The file results/abstracts.txt has been created successfully.\n")
except:
    print("\n An error ocurred while creating the file.\n")
    