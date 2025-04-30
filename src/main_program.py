
import sys, re, os
from peter_fasta_class import Fasta
#from helper_module import get_data_path,load_motif,find_motif
from helper_module import get_data_path,load_motif
from nye_find import check_deviation,find_motif

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
# Load motif and save the function output
motif_list, penalty_list, minimum_gap, maximum_gap = load_motif(motif_file)
print(motif_list, penalty_list, minimum_gap, maximum_gap)

#get matches using find_motif
print("Matches are listed as (start position, penalty score, match)")
print("The header corresponding to the match is printed immediately before the match")
for header, sequence in fasta:
    for match in find_motif(sequence, motif_list, penalty_list, max_deviation, minimum_gap, maximum_gap):
        print(header)
        print(match)


#for match in find_motif(fasta.sequences[1], motif_list, penalty_list, max_deviation):
#        print(match)

