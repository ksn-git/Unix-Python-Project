
import sys, re, os
from peter_fasta_class import Fasta
from helper_module import get_data_path,load_motif

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


def find_motif(sequence, motif_list, penalty_list, max_deviation, minimum_gap, maximum_gap):
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
    for i in range(0, len(sequence) - len(motif_list) - maximum_gap + 1):           # Search until the remaining seq is not long enough to be the motif
        # If no gaps in the motif, window is the same length as the motif
        if minimum_gap == 0 and maximum_gap == 0:
            window = sequence[i:(i + len(motif_list))]                
        # If there are gaps in the motif, window is the motif length + maximum number of gaps in motif. 
        # Subtract 1 because there is a star character in the motif list denoting gaps
        else:
            window = sequence[i:(i + len(motif_list) - 1 + maximum_gap)]                # Window of the same length as the longest possible motif. Len = 29
        # Reset deviation score and search window for motif
        deviation = 0
        for j in range(len(window)):
            # If gap has been reached 
            # Check match for all gap sizes. Yield if below max deviation
            if j == star_index:
                # Range of gaps
                for k in range(minimum_gap, maximum_gap + 1): # range should be the length of the 2nd part of the motif
                    # check match for the length of the 2nd part of the motif
                    for m in range(k, k+len_part_2): # maybe minus 1
                        # If several possible characters
                        if isinstance(motif_list[j+m-k+1], set):
                            print(window[j+m], motif_list[j+m-k+1])    
                            if window[j+m] not in motif_list[j+m-k+1]:  
                                deviation += int(penalty_list[j+m-k+1]) 
                                if deviation > max_deviation:
                                    break
                        # If one possible character 
                        else: 
                            if window[j+m] != motif_list[j+m-k+1]:    
                                deviation += int(penalty_list[j+m-k+1])  
                                # If the deviation is larger than the max, break out of the loop
                                if deviation > max_deviation:
                                    break

            # If gap has been reached, the entire window has already been checked
            elif star_index is not None and j >= star_index:
                # Skip this part if the star/gap has already been encountered
                continue

            elif isinstance(motif_list[j], set):
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


print("Matches are listed as (start position, penalty score, match)")
print("The header corresponding to the match is printed immediately before the match")
for header, sequence in fasta:
    for match in find_motif(sequence, motif_list, penalty_list, max_deviation, minimum_gap, maximum_gap):
        print(header)
        print(match)


#for match in find_motif(fasta.sequences[1], motif_list, penalty_list, max_deviation):
#        print(match)

