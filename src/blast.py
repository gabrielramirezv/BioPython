'''
NAME
    blast

VERSION
    1.0
    
TYPE
    Homework

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/blast.py

DESCRIPTION
    Receives a FASTA file and excute BLAST online to get the 
    alignments with the smallest pvalue.

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
    None
                        
VARIABLES DICTIONARY
    

SEE ALSO
    blastn

'''

# Import libraries
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Get the sequence
sequence = SeqIO.read("data/opuntia1.fasta", format="fasta")

# Execute BLAST
blast_xml = NCBIWWW.qblast("blastn", "nr", sequence.seq)

# Get the information from the xml file
blast_record = NCBIXML.read(blast_xml)

# Get the alignments with pvalue < 0.05 and print them
E_VALUE_THRESH = 0.05
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            print("****Alignment****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
            print("score:",hsp.score)
            print(hsp.query[0:75] + "...")
            print(hsp.match[0:75] + "...")
            print(hsp.sbjct[0:75] + "...")
print("\nDone.")