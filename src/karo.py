
import sys, re, os
from peter_fasta_class import Fasta
from helper_module import get_data_path

# O(1). Constant runtime as the same amount of input is always used. 
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

# Always one sequence file but the amount of entries and their length varies
# Loading and verifying the fasta file 
fasta = Fasta()
# O(n*m), n = sequences, m = lines in the sequences.
fasta.load(fasta_file)   
# O(n), n = sequences.                        
fasta.verify("dnax")     

# Opening the motif description file. Works if the file is in a "data" subfolder. 
# check if file exist is in helper_module 
# O(1)
motif_file = get_data_path(motif)     #this is the fasta file, I've added a reference file     

# O(k), k = lines in the motif file. 
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
            gap_positions = row[1]
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
        elif len(row[0]) > 1:
            multiple_chars = set()
            for char in row[0]:
                multiple_chars.add(char)
            motif_list.append(multiple_chars)
            penalty_list.append(row[1])
        # Lines specifying important positions   
        else:
            motif_list.append(row[0])
            penalty_list.append(row[1])
    infile.close()
    return motif_list, penalty_list, minimum_gap, maximum_gap

# Initialize
motif_list, penalty_list = [], []
minimum_gap, maximum_gap = 0, 0
# Load motif and save the function output
motif_list, penalty_list, minimum_gap, maximum_gap = load_motif(motif_file)
print(motif_list, penalty_list, max_deviation, minimum_gap, maximum_gap)

# Runtime so far
# O(1 + n*m + n + 1 + k) = O(n*m + n + k) = O(n*m + n), n = sequences, m = length of sequences
# The number of lines in the motif file is assumed to be small compared to the number of sequences and their length.

def find_motif(sequence, motif_list, penalty_list, max_deviation, minimum_gap, maximum_gap):
    """ Generator that searches for motif """
    """ Yields position, deviation and sequence when a match is found """
    """ Tries to make a match in each window until max deviation is reached """
    # O(1)
    # Find start position of gap and length of 2nd part of the motif
    try:
        star_index = motif_list.index("*")      
        len_part_2 = len(motif_list) - star_index - 1
    # If no gaps in motif, ignore
    except ValueError:
        star_index = None

    # Start searching for matches
    # O(m), m = sequence length. i will be increase the same amount as the sequence length.
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

        # O(j), j = motif length. j skips the gaps but is dependent on motif length.
        for j in range(len(window)):
            # If gap has been reached 
            # Check match for all gap sizes. Yield if below max deviation
            if j == star_index:
                deviation_from_first_part = deviation
                # Range of gaps
                # O(k), k = gap range. Scales with the number of possible gap sizes.
                for k in range(minimum_gap, maximum_gap + 1): # range should be the length of the 2nd part of the motif
                    # check match for the length of the 2nd part of the motif
                    deviation = deviation_from_first_part
                    if i == 784:
                        print("k",k)
                        print("deviation score", deviation)
                    # O(j/2). Scales with the length of the 2nd part of the motif. I assume the first and second part are equal length
                    for m in range(k, k+len_part_2): # maybe minus 1
                        if i == 784:
                            print("m", m)
                        # If several possible characters
                        if isinstance(motif_list[j+m-k+1], set):
                            if i == 784:
                                print("j+m", str(j+m), "j+m-k+1", str(j+m-k+1))
                                print(window[j+m], motif_list[j+m-k+1])    
                            if window[j+m] not in motif_list[j+m-k+1]:  
                                deviation += int(penalty_list[j+m-k+1]) 
                                if deviation > max_deviation:
                                    break
                        # If one possible character 
                        else: 
                            if i == 784:
                                print("j+m", str(j+m), "j+m-k+1", str(j+m-k+1))
                                print(window[j+m], motif_list[j+m-k+1])
                            if window[j+m] != motif_list[j+m-k+1]:    
                                deviation += int(penalty_list[j+m-k+1])  
                                # If the deviation is larger than the max, break out of the loop
                                if deviation > max_deviation:
                                    break
                        
                        # After testing a gap, yield match if the deviation is below the max
                        # O(k). Scales with the number of possible gap sizes since the check has to be done for each.
                        if m == k + len_part_2 - 1:
                            if deviation <= max_deviation:
                                # Adjust the printed window to only show the matched part of the sequence
                                if k != maximum_gap:
                                    substract_from_window = maximum_gap - k
                                    adjusted_window = window[:-substract_from_window]
                                    yield((i, deviation, window[0:star_index]+"*"+ adjusted_window[-len_part_2:]))
                                else:
                                    yield((i, deviation, window[0:star_index]+"*"+ window[-len_part_2:]))
                                

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
            # O(j). Scales with the length of the motif as the comparison has to be done for each character in the motif.
            else:
                if window[j] != motif_list[j]:
                    deviation += int(penalty_list[j])
                    # If the deviation is larger than the max, break out of the loop
                    if deviation > max_deviation:
                        break

        # Yielding matches for sequences without gaps
        if deviation <= max_deviation and star_index is None:
            yield((i, deviation, window))                           # Return the position, deviation and match


# O(n). Scales with number of sequences you put into find_motif. The actual runtime of find_motif is analyzed above.
print("Matches are listed as (start position, penalty score, match)")
print("The header corresponding to the match is printed immediately before the match")
for header, sequence in fasta:
    for match in find_motif(sequence, motif_list, penalty_list, max_deviation, minimum_gap, maximum_gap):
        print(header)
        print(match)

# Runtime for the entire program

# O(n*m + 1 + m + j + k + j/2 + k + j + n) = O(n*m + m + j + k + n)
# The number of gap sizes, k, is assumed to be very small compared to m, the length of sequences.
# The same assumption is made for j, the length of the motif.
# Final runtime is therefore: 
# O(n*m + m + n)


#for match in find_motif(fasta.sequences[1], motif_list, penalty_list, max_deviation):
#        print(match)

