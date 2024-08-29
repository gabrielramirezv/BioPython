'''
NAME
    blast

VERSION
    0.1
    
TYPE
    Homework

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/blast.py

DESCRIPTION
    Receives a FASTA file and excute BLAST online to get the 
    alignements with the smallest pvalue.

CATEGORY
    BLAST
    
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
    python src/blast.py

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
    blastn

'''

# Import libraries

# Execute BLAST

# Get the alignments with pvalue < 0.05 and print them
