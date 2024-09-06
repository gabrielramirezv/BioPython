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
        |- src       # The script must be executed from this directory
USAGE
    python src/blast.py

ARGUMENTS
    None
                        
VARIABLES DICTIONARY
    
SEE ALSO
    blastn
'''

# Import libraries
import argparse
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Create the parser
parser = argparse.ArgumentParser(description="Receives a FASTA file "
                                 + "and execute BLAST online to get " 
                                 + "the alignments with the smallest " 
                                 + "pvalue.")

# Define the arguments
parser.add_argument("-i", "--inputfile",
                    help="Path to the file with the sequence to blast",
                    type=str,
                    required=True)

parser.add_argument("-f", "--format",
                    help="Format of the sequence file",
                    type=str,
                    required=False,
                    default="fasta")

parser.add_argument("-b", "--blasttype",
                    help="Blast type to be executed",
                    type=str,
                    required=False,
                    default="blastn")

parser.add_argument("-d", "--database",
                    help="Database to blast",
                    type=str,
                    required=False,
                    default="nr")

parser.add_argument("-e", "--evalue",
                    help="e value limit for the alignments",
                    type=float,
                    required=False,
                    default=0.05)

# Get the parameters to run BLAST
args = parser.parse_args()
seq_file = args.inputfile
input_format = args.format
blast_type = args.blasttype
database = args.database
E_VALUE_THRESH = args.evalue

# Read the sequences from the input file
for seq in SeqIO.parse(seq_file, format=input_format):
    # Execute BLAST
    blast_xml = NCBIWWW.qblast(blast_type, database, seq.seq)

    # Parse the BLAST XML
    blast_record = NCBIXML.read(blast_xml)

    # Get the best alignments with pvalue < 0.05 and print them
    E_VALUE_THRESH = 0.05
    for alignment in blast_record.alignments:
        hsp_count = 0
        if alignment.hsps[hsp_count].expect < E_VALUE_THRESH:
            print(f"****Best Alignment for {seq.id}****")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", alignment.hsps[hsp_count].expect)
            print("score:",alignment.hsps[hsp_count].score)
            print(alignment.hsps[hsp_count].query[0:75] + "...")
            print(alignment.hsps[hsp_count].match[0:75] + "...")
            print(alignment.hsps[hsp_count].sbjct[0:75] + "...\n")
            break
        else:
            hsp_count += 1
print("\nDone.\n")
