'''
NAME
    lineage_comparison

VERSION
    1.0

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/lineage_comparison.py

DESCRIPTION
    Search taxonomy files for two different species and compare their 
    lineages

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
    python src/lineage_comparison.py

ARGUMENTS
    None

SEE ALSO
    Esearch

'''

# Import libraries
import argparse
from Bio import Entrez

# Create parser
parser = argparse.ArgumentParser(description="Search taxonomy files " 
                                 + "for two different species and " 
                                 + "compare their lineages")

# Store the arguments
parser.add_argument("-e", "--email",
                    help="Email to use in Entrez",
                    type=str,
                    required=False,
                    default="gramirez@lcg.unam.mx")

parser.add_argument("-o", "--organism1",
                    help="First organism to search",
                    type=str,
                    required=False,
                    default="Notoryctes typhlops")

parser.add_argument("-p", "--organism2",
                    help="Second organism to search",
                    type=str,
                    required=False,
                    default="Chrysochloris asiatica")

args = parser.parse_args()
email = args.email
organism_1 = args.organism1
organism_2 = args.organism2


# Register an email account
Entrez.email = email

# Get the ID for Notoryctes typhlops
with Entrez.esearch(db="Taxonomy", 
                    term=organism_1) as handle:
    record = Entrez.read(handle)
    organism_taxonomy_id = record["IdList"]

# Get taxonomy file for Notoryctes typhlops
with Entrez.efetch(db = "Taxonomy", 
                   id = organism_taxonomy_id, 
                   retmode="xml") as handle:
    organism_1_taxonomy = Entrez.read(handle)

# Search the information for Chrysochloris asiatica
with Entrez.esearch(db="Taxonomy", 
                    term=organism_2) as handle:
    record = Entrez.read(handle)
    organism_taxonomy_id = record["IdList"]

# Get taxonomy file for Chrysochloris asiatica
with Entrez.efetch(db = "Taxonomy", 
                   id = organism_taxonomy_id, 
                   retmode="xml") as handle:
    organism_2_taxonomy = Entrez.read(handle)

# Get lineages of both organisms
print(f"{organism_1} lineage: " 
      + "{organism_1_taxonomy[0]['Lineage'].replace('; ', ', ')}")
print(f"\n{organism_2} lineage: {organism_2_taxonomy[0]['Lineage']}")

# Compare lineages
organism_1_tax_levels = organism_1_taxonomy[0]["Lineage"].split("; ")
organism_2_tax_levels = organism_2_taxonomy[0]["Lineage"].split("; ")

# Report results to user
different_elements = list()
print("\nThese are the common elements in the lineages of both organisms:")
if len(organism_1_tax_levels) > len(organism_2_tax_levels):
    for level in organism_1_tax_levels:
        if level in organism_2_tax_levels:
            print(level)
        else:
            different_elements.append(level)
    print("\nThese are the different elements between both organisms " 
          + "lineages:")
    for element in different_elements:
        print(element)
else:
    for level in organism_2_tax_levels:
        if level in organism_1_tax_levels:
            print(level)
        else:
            different_elements.append(level)
    print("\nThese are the different elements between both organisms:")
    for element in different_elements:
        try:
            print(f"{element}\t" 
                + organism_1_tax_levels[organism_2_tax_levels.index(element)])
        except:
            print(f"{element}")
