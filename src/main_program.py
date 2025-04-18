
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
motif_file = get_data_path('motif.fsa')     #this is the fasta file, I've added a reference file     
infile = open(motif_file,'r')

# Arranging motif into two lists
# Motif list:   ATCGGATC******AGTCGTTA
# Penalty list: 8767894200000087698432
motif_list = []
penalty_list = []
for line in infile:
    # Skip descriptive lines
    if line.startswith("#"):
        continue
    else:
        # Line specifying unimportant positions
        if line.startswith("*"):
            unimportant_positions = line.strip().split(sep="\t")[1]
            # If there is a defined number of unimportant positions
            try:
                unimportant_positions = int(unimportant_positions)
                for i in range(0, unimportant_positions):
                    motif_list.append("*")
                    penalty_list.append(0)
                    continue 
            except ValueError:
                pass                
            # If there is a range of unimportant positions
            if isinstance(unimportant_positions, str):
                result = re.search(r"(\d+)-(\d+)", unimportant_positions)
                if result is not None:
                    minimum = int(result.group(1))
                    maximum = int(result.group(2))
                    for i in range(0, minimum):
                        motif_list.append("*")
                        penalty_list.append(0)
                else:
                    raise ValueError("Invalid unimportant positions format")
                    sys.exit(1)
        # Lines specifying important positions
        else:
            motif_list.append(line.strip().split(sep="\t")[0])
            penalty_list.append(line.strip().split(sep="\t")[1])
infile.close()

#find motif
#I think this should be added to the class...
def find_motif(sequence,motif_list,penalty_list,max_deviation):
    """Generator that yields a matching motif."""
    #to avoid index error
    print('placeholder for now')

    ##psudocode
    
    # for seq in sequences
        #reset mismatch score
        #reset temp match 

        # for each sequence position in seq (using sliding window)
            #check if sequence position value is below deviation
            # if temp match != '' and last postion in match is reached (x positions filled with not '*'?)
                #save/yield? (temp match,position,mismatch score) to match list for the first entry
                #reset mismatch score
                #reset temp match
            # if over deviation,
                #move to next start position
                #reset mismatch score
                #reset temp match
            # if below deviation,
                #add possible mismatch value to mismatch score
                #add position to temp match
            #if match is in progress and skip position is reached (start to 15)
                #add '*' to file and don't 
    
# need to find out how to jump back to last starting position. 
# might be missing some conditions


 
def find_motif(sequence, motif_list, penalty_list, max_deviation):
    """ Generator that searches for motif """
    """ Yields position, deviation and sequence when a match is found """
    """ Tries to make a match in each window until max deviation is reached """

    # OBS: ONLY WORKS WHEN EACH POSITION HAS ONE POSSIBLE LETTER
    for i in range(0, len(sequence) - len(motif_list)):             # Search until the remaining seq is not long enough to be the motif
        window = sequence[i:i + len(motif_list)]                    # Window of the same length as the motif
        print(window)
        deviation = 0
        for j in range(len(window)):
            # If the window is not equal to the motif, add penalty score
            if window[j] != motif_list[j]:
                deviation += int(penalty_list[j])
                # If the deviation is larger than the max, break out of the loop
                if deviation > max_deviation:
                    break
        # If the deviation is less than the max, print the match
        if deviation <= max_deviation:
            yield((i, deviation, window))         # Return the position, deviation and match

find_motif(fasta.sequences[0], motif_list, penalty_list, max_deviation)

# Searching for the motif in each entry in the fasta file
#for sequence in fasta.sequences:
#    find_motif(motif_list, penalty_list, sequence, max_deviation)



