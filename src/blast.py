'''
NAME
    blast

VERSION
    3.0
    
TYPE
    Homework

AUTHOR
    Gabriel Ramirez Vilchis

GITHUB
    https://github.com/gabrielramirezv/BioPython/blob/master/src/blast.py

DESCRIPTION
    Receives a FASTA file and excute BLAST online to get the 
    alignements with the smallest e value. If the file has only one 
    sequence, the program will print the alignments with an e value 
    lower than 0.05 or the value entered by the user. If the file 
    contains many sequences, the program will print only the best hit
    according to the score for each sequence.

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
    python blast.py [-h] -i INPUTFILE [-f FORMAT] [-b BLASTTYPE] 
        [-d DATABASE] [-e EVALUE] [-x]

ARGUMENTS
    -h, --help           show this help message and exit
    -i INPUTFILE, --inputfile INPUTFILE
                        Path to the file with the sequence to blast
    -f FORMAT, --format FORMAT
                        Format of the sequence file
    -b BLASTTYPE, --blasttype BLASTTYPE
                        Blast type to be executed
    -d DATABASE, --database DATABASE
                        Database to blast
    -e EVALUE, --evalue EVALUE
                        e value limit for the alignments
    -x, --excludeitself
                        If selected, the same sequence will not be 
                        able to be the best hit for itself, when 
                        analizing a file with many sequences
    
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
                                 + "e value.")

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

parser.add_argument("-x", "--excludeitself",
                    help="If selected, the same sequence will not be "
                    + "able to be the best hit for itself, when "
                    + "analizing a file with many sequences",
                    action="store_true",
                    default=False)

# Get the parameters to run BLAST
args = parser.parse_args()
seq_file = args.inputfile
input_format = args.format
blast_type = args.blasttype
database = args.database
exclude_itself = args.excludeitself
E_VALUE_THRESH = args.evalue

# Try to parse the sequences from the file
try:
    # Read the sequences from the input file
    seqs = list(SeqIO.parse(seq_file, format=input_format))

# If the file does not exist, ask the user to try again
except:
    print("The file was not found, please try to write the path of "
          + "the input file again.")

else:
    # This program just prints the best hit for each sequence if there 
    # is more than one sequence to analyze
    if (len(seqs) > 1) and not exclude_itself:

        # Search the best hit for each sequence 
        for seq in seqs:

            # Execute BLAST
            blast_xml = NCBIWWW.qblast(blast_type, database, seq.seq)

            # Read the BLAST XML
            blast_record = NCBIXML.read(blast_xml)

            # Get the best alignments with a lower e value than the 
            # thresh and print them
            for alignment in blast_record.alignments:
                # Since the alignments are sorted by score, iterate 
                # the alignments until the program finds an e value
                # lower than the thresh
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

    elif (len(seqs) > 1) and exclude_itself:
        
        # Search the best hit for each sequence 
        for seq in seqs:

            # Execute BLAST
            blast_xml = NCBIWWW.qblast(blast_type, database, seq.seq)

            # Read the BLAST XML
            blast_record = NCBIXML.read(blast_xml)

            # Get the best alignments with a lower e value than the thresh 
            # and print them
            for alignment in blast_record.alignments:
                # Since the alignments are sorted by score, iterate the 
                # alignments until the program finds an e value lower than 
                # the thresh
                hsp_count = 0
                if ((alignment.hsps[hsp_count].expect < E_VALUE_THRESH) 
                    and (seq.id not in alignment.title)):
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

    # If there is just a single sequence to analyze, the program will 
    # print the alignments with a lower e value than the thresh
    else:

        # While there is an only sequence, it will be stored in a 
        # variable for legibility
        seq = seqs[0]

        # Execute BLAST
        blast_xml = NCBIWWW.qblast(blast_type, database, seq.seq)

        # Read the BLAST XML
        blast_record = NCBIXML.read(blast_xml)

        # Evaluate the e evalues and print the alignments under the 
        # e value limit
        number_of_alignment = 0 
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    number_of_alignment += 1
                    print(f"****Alignment {number_of_alignment}****")
                    print("sequence:", alignment.title)
                    print("length:", alignment.length)
                    print("e value:", hsp.expect)
                    print("score:",hsp.score)
                    print(hsp.query[0:75] + "...")
                    print(hsp.match[0:75] + "...")
                    print(hsp.sbjct[0:75] + "...\n")

    # Let the user know that the process is finished
    print("\nDone.\n\n")
