
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
motif_file = get_data_path(motif)          

infile = open(motif_file,'r')
# Arranging motif into two lists
# Motif list:   ATCGGATC*AGTCGTTA
# Penalty list: 87678942087698432
motif_list = []
penalty_list = []
for line in infile:
    # Skip descriptive lines
    if line.startswith("#"):
        continue
    # Line specifying unimportant positions
    elif line.startswith("*"):
        unimportant_positions = line.strip().split(sep="\t")[1]
        # If there is a defined number of unimportant positions
        try:
            unimportant_positions = int(unimportant_positions)
            minimum = unimportant_positions
            maximum = unimportant_positions
            motif_list.append("*")
            penalty_list.append(0) 
        except ValueError:
            pass                
        # If there is a range of unimportant positions
        if isinstance(unimportant_positions, str):
            result = re.search(r"(\d+)-(\d+)", unimportant_positions)
            if result is not None:
                minimum = int(result.group(1))
                maximum = int(result.group(2))
                motif_list.append("*")
                penalty_list.append(0)
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

"""
#saved use case code for find_motif
print("Matches are listed as (start position, penalty score, match)")
print("The header corresponding to the match is printed immediately before the match")
for header, sequence in fasta:
    for match in find_motif(sequence, motif_list, penalty_list, max_deviation):
        print(header)
        print(match)

#for match in find_motif(fasta.sequences[1], motif_list, penalty_list, max_deviation):
#        print(match)
"""
