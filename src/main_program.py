
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
print(motif_list,penalty_list)
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
    # Find start position of gap and length of 2nd part of the motif
    try:
        star_index = motif_list.index("*")      
        len_part_2 = len(motif_list) - star_index - 1
    # If no gaps in motif, ignore
    except ValueError:
        star_index = None

    # Start searching for matches
    for i in range(0, len(sequence) - len(motif_list)):             # Search until the remaining seq is not long enough to be the motif
        window = sequence[i:(i + len(motif_list) + maximum - 1)]    # Window of the same length as the longest possible motif. Len = 29
        deviation = 0
        for j in range(len(window)):
              
            # If gap has been reached 
            # check deviation for 15 gaps, 16 gaps, 17 gaps. Save deviation
            # Choose gaps corresponding to lowest deviation
            # Continue searching for match after chosen gap size

            # skip min amount of starts, then min +1 ... up to max amount of starts
            # motif list should have one star 
            # while motif life is still long enough, check if equal to window

            # OBS: Motif is currently 29. But j + m(max) = 29, so index error
            # For gap = 17, k = 17, m = 17-23, j = 6. these are human indices/lengths, not python indices! Change this in code
            # If the motif contains a gap / star and it has been reached
            if j == star_index:
                # Range of gaps 
                for k in range(minimum, maximum + 1): # range should be the length of the 2nd part of the motif
                    # check match for the length of the 2nd part of the motif
                    for m in range(k, k+len_part_2):
                        # If several possible characters
                        if isinstance(motif_list[j+m], set):
                            if window[j+m] not in motif_list[j+m]:
                                deviation += int(penalty_list[j])
                                if deviation > max_deviation:
                                    break
                        # If one possible character
                        else: 
                            if window[j+m] != motif_list[j+m]:
                                deviation += int(penalty_list[j])
                                # If the deviation is larger than the max, break out of the loop
                                if deviation > max_deviation:
                                    break
                    if deviation <= max_deviation:
                        yield((i, deviation, window))

            # If several possible characters
            if isinstance(motif_list[j], set):
                # If the character is not in the list of possible characters, add penalty score
                if window[j] not in motif_list[j]:
                    deviation += int(penalty_list[j])
                    # If the deviation is larger than the max, break out of the loop
                    if deviation > max_deviation:
                        break

            # If one possible character
            else:
                if window[j] != motif_list[j]:
                    deviation += int(penalty_list[j])
                    # If the deviation is larger than the max, break out of the loop
                    if deviation > max_deviation:
                        break

        # If the deviation is less than the max, yield the match
        if deviation <= max_deviation:
            yield((i, deviation, window))                           # Return the position, deviation and match

# for sequence in fasta.sequences:
for match in find_motif(fasta.sequences[-5], motif_list, penalty_list, max_deviation):
    print(match)


