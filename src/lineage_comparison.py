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
                        
VARIABLES DICTIONARY
    

SEE ALSO
    Esearch

'''

# Import libraries
from Bio import Entrez

# Register an email account
Entrez.email = "gramirez@lcg.unam.mx"

# Get the ID for Notoryctes typhlops
with Entrez.esearch(db="Taxonomy", 
                    term="Notoryctes typhlops") as handle:
    record = Entrez.read(handle)
    organism_taxonomy_id = record["IdList"]

# Get taxonomy file for Notoryctes typhlops
with Entrez.efetch(db = "Taxonomy", 
                   id = organism_taxonomy_id, 
                   retmode="xml") as handle:
    NotoryctesTyphlops_taxonomy = Entrez.read(handle)

# Search the information for Chrysochloris asiatica
with Entrez.esearch(db="Taxonomy", 
                    term="Chrysochloris asiatica") as handle:
    record = Entrez.read(handle)
    organism_taxonomy_id = record["IdList"]

# Get taxonomy file for Chrysochloris asiatica
with Entrez.efetch(db = "Taxonomy", 
                   id = organism_taxonomy_id, 
                   retmode="xml") as handle:
    ChrysochlorisAsiatica_taxonomy = Entrez.read(handle)

# Get lineages of both organisms
print(f"Notoryctes typhlops lineage: {NotoryctesTyphlops_taxonomy[0]['Lineage'].replace('; ', ', ')}")
print(f"\nChrysochloris asiatica lineage: {ChrysochlorisAsiatica_taxonomy[0]['Lineage']}")

# Compare lineages
NotoryctesTyphlops_tax_levels = NotoryctesTyphlops_taxonomy[0]["Lineage"].split("; ")
ChrysochlorisAsiatica_tax_levels = ChrysochlorisAsiatica_taxonomy[0]["Lineage"].split("; ")

# Report results to user
different_elements = list()
print("\nThese are the the common elements in the lineages of both organisms:")
if len(NotoryctesTyphlops_tax_levels) > len(ChrysochlorisAsiatica_tax_levels):
    for level in NotoryctesTyphlops_tax_levels:
        if level in ChrysochlorisAsiatica_tax_levels:
            print(level)
        else:
            different_elements.append(level)
    print("\nThese are the different elements between both organisms lineages:")
    for element in different_elements:
        print(element)
else:
    for level in ChrysochlorisAsiatica_tax_levels:
        if level in NotoryctesTyphlops_tax_levels:
            print(level)
        else:
            different_elements.append(level)
    print("\nThese are the different elements between both organisms:")
    for element in different_elements:
        try:
            print(f"{element}\t{NotoryctesTyphlops_tax_levels[ChrysochlorisAsiatica_tax_levels.index(element)]}")
        except:
            print(f"{element}")
