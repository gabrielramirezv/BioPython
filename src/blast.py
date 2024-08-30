'''
NAME
    blast

VERSION
    2.0
    
TYPE
    Homework

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/blast.py

DESCRIPTION
    Receives a FASTA file and excute BLAST online to get the 
    alignments with the smallest evalue.

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
    From the workspace main place with the structure described above 
    in front of you:
        python src/blast.py [-h] -i INPUTFILE [-f FORMAT] 
            [-b BLASTTYPE] [-d DATABASE] [-o OUTPUTFILE] 
            [-e EVALUE]

ARGUMENTS
    -h, --help              Show this help message and exit
    -i INPUTFILE, --inputfile INPUTFILE
                            Path to the file with the sequence to blast
    -f FORMAT, --format FORMAT
                            Format of the sequence file
    -b BLASTTYPE, --blasttype BLASTTYPE
                            Blast type to be executed
    -d DATABASE, --database DATABASE
                            Database to blast
    -o OUTPUTFILE, --outputfile OUTPUTFILE
                            Path to the file to create with the results
    -e EVALUE, --evalue EVALUE
                            e value limit for the alignments
                        
VARIABLES DICTIONARY
    args: Arguments received from the user
    blast_type: Type of blast to execute
    blast_xml: XML result of the BLAST
    blst_record: Content of the BLAST XML after being parsed
    data_base: Data base of BLAST to be aligned by the query sequence
    E_VALUE_THRESH: e value limit for the alignments
    file_format: Format of the input file
    number_of_alignment: Count of the alignments being reported
    output_file: Path to the file to create with the results
    parser: Arguments parser
    sequence: Sequence after being read from the input file
    sequence_file: Path to the file with the sequence to blast
    
SEE ALSO
    blastn

'''

# Import libraries
import argparse
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

# Define function to print the results or write them in a file
def print_results(number_of_alignment, alignment, e_value):
    '''
    It writes each alignment BLAST results in a document or print them
    to sreen, according to the user's choice.
        Parameters:
            number of alignment (int): Number of the alignment to be 
                written
            alignment (obj): Alignment to be reported
            hsp.expect (float): e value of the alignment
    '''
    if output_file:
        with open(output_file, 'w') as output:
            output.write(f"***Alignment {number_of_alignment}***"
                         + f"\nsequence: {alignment.title}"
                         + f"\nlength: {alignment.length}"
                         + f"\ne value: {e_value}")
    else:
        print(f"***Alignment {number_of_alignment}***\n"
              "\nsequence:" , alignment.title,
              "\nlength: ", alignment.length,
              "\ne value: ", e_value)

# Create parser
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

parser.add_argument("-o", "--outputfile",
                    help="Path to the file to create with the results",
                    type=str,
                    required=False)

parser.add_argument("-e", "--evalue",
                    help="e value limit for the alignments",
                    type=float,
                    required=False,
                    default=0.05)

# Store the arguments
args = parser.parse_args()
sequence_file = args.inputfile
file_format  = args.format
output_file = args.outputfile
blast_type = args.blasttype
data_base = args.database
E_VALUE_THRESH = args.evalue

# Read the sequence from the input file
sequence = SeqIO.read(sequence_file, format=file_format)

# Execute BLAST
blast_xml = NCBIWWW.qblast(blast_type, data_base, sequence.seq)

# Get the information from the xml file
blast_record = NCBIXML.read(blast_xml)

# Get the alignments with evalue < 0.05 and print them
number_of_alignment = 0
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < E_VALUE_THRESH:
            number_of_alignment += 1
            print_results(number_of_alignment, alignment, hsp.expect)
            
# Let the user know that the process is finished
print("\nDone.")
