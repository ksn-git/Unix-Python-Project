
import sys, re, os
from peter_fasta_class import Fasta
from helper_module import get_data_path

# Getting filename with motif from command line
if len(sys.argv) != 4:                              
    if len(sys.argv) == 1:
        motif = input("Please enter the name of the file with the motif to look for: ")
        fasta_file = input("Please enter the name of the fasta file to search in: ")
        max_deviation = int(input("Please enter the maximum deviation (int): "))
    else:
        print("You must supply a file name with a motif to look for")
        print("Usage: <main_program.py> <motif> <fastafile> <max deviation>")
        sys.exit(1)
else:
    motif = sys.argv[1]                              
    fasta_file = sys.argv[2]
    max_deviation = int(sys.argv[3])

# Loading and verifying the fasta file 
fasta = Fasta()
fasta.load(fasta_file)                           
fasta.verify("dnax")     

# Opening the motif description file. Works if the file is in a "data" subfolder. 
# check if file exist is in helper_module 
motif_file = get_data_path(motif)     #this is the fasta file, I've added a reference file     

def load_motif(motif_file):
    """ Returns a list of the motif and associated penalty scores as well as minimum and maximum number of gaps """
    # Arranging motif into two lists. Example output:
    # Motif list:   ATCGGATC*AGTCGTTA
    # Penalty list: 87678942047698432
    
    # Open file with motif
    infile = open(motif_file,'r')
    # Initialize objects
    motif_list, penalty_list = [], []
    minimum_gap, maximum_gap = 0, 0
    
    for line in infile:
        row = line.strip().split(sep="\t")
        # Skip descriptive lines
        if line.startswith("#"):
            continue
        # Line specifying unimportant positions / gap in motif
        elif line.startswith("*"):
            gap_positions = line.strip().split(sep="\t")[1]
            # If there is a defined number of unimportant positions
            try:
                gap_positions = int(gap_positions)
                minimum_gap = gap_positions
                maximum_gap = gap_positions
                motif_list.append("*")
                penalty_list.append(0) 
            except ValueError:       
                # If there is a range of unimportant positions
                    result = re.search(r"(\d+)-(\d+)", gap_positions)
                    if result is not None:
                        minimum_gap = int(result.group(1))
                        maximum_gap = int(result.group(2))
                        motif_list.append("*")
                        penalty_list.append(0)
                    # If not an integer or a range, then the format is invalid
                    else:
                        raise ValueError("Invalid unimportant positions format")
                        sys.exit(1)
        # Lines with >1 possible character 
        elif len(line.strip().split(sep="\t")[0]) > 1:
            multiple_chars = set()
            for char in line.strip().split(sep="\t")[0]:
                multiple_chars.add(char)
            motif_list.append(multiple_chars)
            penalty_list.append(line.strip().split(sep="\t")[1])
        # Lines specifying important positions   
        else:
            motif_list.append(line.strip().split(sep="\t")[0])
            penalty_list.append(line.strip().split(sep="\t")[1])
    infile.close()


"""print("Matches are listed as (start position, penalty score, match)")
print("The header corresponding to the match is printed immediately before the match")
for header, sequence in fasta:
    for match in find_motif(sequence, motif_list, penalty_list, max_deviation):
        print(header)
        print(match)

#for match in find_motif(fasta.sequences[1], motif_list, penalty_list, max_deviation):
#        print(match)
"""
