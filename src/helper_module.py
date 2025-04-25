#!/bin/env python3

#dependencies add_to_sys_path and get_data_path
import os
import sys
#dependency detect_motif
import re

#change sys_path
def add_to_sys_path(relative_path_from_root: str):
    """Add directory to sys path from project root."""

    #check if relative_path_from_root is empty string 
    if not relative_path_from_root:
        raise ValueError(f'The relative path is empty: {relative_path_from_root}')

    #find current path
    current_file = os.path.abspath(__file__)
    #go back and look at project root
    project_root = os.path.dirname(os.path.dirname(current_file))
    #get target path
    target_path = os.path.join(project_root,relative_path_from_root)

    #check if directory exist
    if not os.path.isdir(target_path):
        raise FileNotFoundError(f'The directory {target_path} does not exist.')
    #add target_path if needed
    if target_path not in sys.path:
        sys.path.insert(0,target_path)
    

#add file path ('data' is the default directive for the data)
def get_data_path(filename,data = 'data'):
    """Finds full data path in directory"""

    #find current path
    current_file = os.path.abspath(__file__)
    #go back and look at project root
    project_root = os.path.dirname(os.path.dirname(current_file))
    #return full path
    full_path = os.path.join(project_root, data, filename)

    #check if file exist
    if not os.path.isfile(full_path):
        raise FileNotFoundError(f'File {filename} does not exist at path {full_path}')
    
    return full_path

#load motif
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
    return motif_list, penalty_list, minimum_gap, maximum_gap



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

            # skip min amount of stars, then min +1 ... up to max amount of starts
            # motif list should have one star 
            # while motif list is still long enough, check if equal to window

            # For gap = 17: k = 17, m = [17-23], j = 6. these are human indices/lengths, not python indices! Change this in code
            # If the motif contains a gap / star and it has been reached
            # If TTCAGA*, then j = 6 when star is reached
            if j == star_index:
                # Range of gaps, e.g. 15, 16, 17
                for k in range(minimum_gap, maximum_gap + 1): # range should be the length of the 2nd part of the motif
                    # check match for the length of the 2nd part of the motif
                    for m in range(k, k+len_part_2 - 1):
                        print(f"Checking window index: {j + m}")
                        # If several possible characters
                        if isinstance(motif_list[j+m-k+1], set):    
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
            yield((i, deviation, window))