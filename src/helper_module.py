#!/bin/env python3

#dependencies add_to_sys_path and get_data_path
import os
import sys
#dependency load_motif
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
                        pos1 = int(result.group(1))
                        pos2 = int(result.group(2))
                        minimum_gap = min(pos1,pos2)
                        maximum_gap = max(pos1,pos2)
                        if pos1 != minimum_gap or pos2 != maximum_gap:
                            print(f'The position of the minimum and maximum gap has been switched from [{pos1}:{pos2}] to [{minimum_gap}:{maximum_gap}]')

                        motif_list.append("*")
                        penalty_list.append(0)
                    # If not an integer or a range, then the format is invalid
                    else:
                        raise ValueError("Invalid unimportant positions format")
        # Lines with >1 possible character 
        elif len(row[0]) > 1:
            multiple_chars = set()
            for char in row[0]:
                multiple_chars.add(char)
            motif_list.append(multiple_chars)
            try: penalty_list.append(int(row[1]))
            except ValueError as err: 
                raise ValueError(f"Penalty must be an integer, got '{row[1]}' instead.") from err                
        # Lines specifying important positions   
        else:
            motif_list.append(row[0])
            try: penalty_list.append(int(row[1]))
            except ValueError as err: 
                raise ValueError(f"Penalty must be an integer, got '{row[1]}' instead.") from err                
    infile.close()
    return motif_list, penalty_list, minimum_gap, maximum_gap

def find_motif(sequence, motif_list, penalty_list, max_deviation, minimum_gap, maximum_gap):
    """ Generator that searches for motif """
    """ Yields position, deviation and sequence when a match is found """
    """ Tries to make a match in each window until max deviation is reached """
    # Check input
    if len(motif_list) != len(penalty_list) or motif_list == [] or penalty_list == []:
        raise ValueError('The list of the motif or the penalty is not the same length or both are empty.')

    # Find start position of gap and length of 2nd part of the motif
    try:
        star_index = motif_list.index("*")      
        len_part_2 = len(motif_list) - star_index - 1
    # If no gaps in motif, ignore
    except ValueError:
        star_index = None

    #calculate window_size
    if minimum_gap == 0 and maximum_gap == 0:
        window_size = len(motif_list) 
    else:
        window_size = len(motif_list) - 1 + maximum_gap

    # Start searching for matches
    for i in range(0, len(sequence) - window_size + 1):           # Search until the remaining seq is not long enough to be the motif
        # Reset deviation score and search window for motif
        window = sequence[i:i + window_size]
        deviation = 0
        match_str = ''

        for j in range(len(window)):
            # If gap has been reached 
            # Check match for all gap sizes. Yield if below max deviation
            if j == star_index:
                deviation_from_first_part = deviation
                match_str_1 = match_str
                # Range of gaps
                for k in range(minimum_gap, maximum_gap + 1): # range should be the length of the 2nd part of the motif
                    # check match for the length of the 2nd part of the motif
                    deviation = deviation_from_first_part
                    match_str_2 = ''
                    for m in range(k, k+len_part_2): # maybe minus 1
                        motif_index = j + m - k + 1
                        if m >= len(window):
                            break       #avoids index out of range
                        # If several possible characters
                        if isinstance(motif_list[motif_index], set):    
                            if window[j+m] not in motif_list[motif_index]:  
                                deviation += int(penalty_list[motif_index]) 
                                match_str_2 += 'X'
                                if deviation > max_deviation:
                                    break
                            else:
                                match_str_2 += window[j+m]
                        # If one possible character 
                        else: 
                            if window[j+m] != motif_list[motif_index]:    
                                deviation += int(penalty_list[motif_index]) 
                                match_str_2 += 'X' 
                                # If the deviation is larger than the max, break out of the loop
                                if deviation > max_deviation:
                                    break
                            match_str_2 += window[j+m]
                        
                        # After testing a gap, yield match if the deviation is below the max
                        if m == k + len_part_2 - 1 and deviation <= max_deviation:
                            yield((i, deviation, match_str_1 + "*" + match_str_2))           

            # If gap has been reached, the entire window has already been checked
            elif star_index is not None and j >= star_index:
                # Skip this part if the star/gap has already been encountered
                continue

            elif j < len(motif_list):
                if isinstance(motif_list[j], set):
                    # If the character is not in the list of possible characters, add penalty score
                    if window[j] not in motif_list[j]:
                        deviation += int(penalty_list[j])
                        match_str += 'X'
                        # If the deviation is larger than the max, break out of the loop
                        if deviation > max_deviation:
                            break
                    else:
                        match_str += window[j]
                # If one possible character
                else:
                    if window[j] != motif_list[j]:
                        deviation += int(penalty_list[j])
                        match_str += 'X'
                        # If the deviation is larger than the max, break out of the loop
                        if deviation > max_deviation:
                            break
                    else:
                        match_str += window[j]

        # Yielding matches for sequences without gaps
        if deviation <= max_deviation and star_index is None:
            yield((i, deviation, match_str))                           # Return the position, deviation and match

'''
def find_motif(sequence, motif_list, penalty_list, max_deviation, minimum_gap, maximum_gap):
    """ Generator that searches for motif """
    """ Yields position, deviation and sequence when a match is found """
    """ Tries to make a match in each window until max deviation is reached """
    # Check input
    if len(motif_list) != len(penalty_list) or motif_list == [] or penalty_list == []:
        raise ValueError('The list of the motif or the penalty is not the same length or both are empty.')

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
                deviation_from_first_part = deviation
                # Range of gaps
                for k in range(minimum_gap, maximum_gap + 1): # range should be the length of the 2nd part of the motif
                    # check match for the length of the 2nd part of the motif
                    deviation = deviation_from_first_part
                    for m in range(k, k+len_part_2): # maybe minus 1
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
                        
                        # After testing a gap, yield match if the deviation is below the max
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
            else:
                if window[j] != motif_list[j]:
                    deviation += int(penalty_list[j])
                    # If the deviation is larger than the max, break out of the loop
                    if deviation > max_deviation:
                        break

        # Yielding matches for sequences without gaps
        if deviation <= max_deviation and star_index is None:
            yield((i, deviation, window))                           # Return the position, deviation and match

'''