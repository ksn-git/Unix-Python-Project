#!/bin/env python3

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