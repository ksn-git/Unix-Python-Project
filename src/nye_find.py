#!/bin/env python3
'''
#check deviation
def check_deviation(window_elem,motif_elem,penalty_elem,deviation,max_deviation):
    """Compare the element from the window with the motif element to return the updated deviation
    and boolean indicating success or not (exceeding max)"""

    #if multiple charactors are found
    if isinstance(motif_elem,set):
        if window_elem not in motif_elem:
            deviation += penalty_elem
            if deviation > max_deviation:
                return deviation, True,'X'      #signal stop now
            return deviation, False, window_elem
    #if a single charactor is found
    else: 
        if window_elem != motif_elem:
            deviation += penalty_elem
            if deviation > max_deviation:
                return deviation, True, 'X'
            return deviation, False, window_elem
    #return False if match continues
    return deviation, False, window_elem         #continue match

def find_motif(sequence, motif_list, penalty_list, max_deviation, minimum_gap, maximum_gap):
    """Generator that searches for motif.
    Yields (position, deviation, matched_sequence) when a match is found."""

    #check if motif and penalty length match
    if len(motif_list) != len(penalty_list) or motif_list == []:
        raise ValueError('The list of the motif or the penalty is not the same length or both are empty.')

    try:
        star_index = motif_list.index("*")
        len_part_2 = len(motif_list) - star_index - 1
    except ValueError:
        star_index = None

    #define motif length (gap or no gap)
    motif_core_length = len(motif_list) - (1 if star_index is not None else 0)

    for i in range(0, len(sequence) - motif_core_length - maximum_gap + 1):
        deviation = 0
        success = False

        if star_index is None:
            # No gaps in motif
            window = sequence[i:i + len(motif_list)]
            match = []

            for j in range(len(motif_list)):
                w = window[j]
                m = motif_list[j]
                p = penalty_list[j]
                deviation, exceeded, char = check_deviation(w, m, p, deviation, max_deviation)
                match.append(char)
                if exceeded:
                    break
            else:
                success = True
                matched_seq = ''.join(match)

        else:
            # Gapped motif
            for gap_size in range(minimum_gap, maximum_gap + 1):
                window = sequence[i:i + len(motif_list) - 1 + gap_size]
                match = []
                deviation = 0
                exceeded = False

                # Before the gap
                for j in range(star_index):
                    w = window[j]
                    m = motif_list[j]
                    p = penalty_list[j]
                    deviation, exceeded, char = check_deviation(w, m, p, deviation, max_deviation)
                    match.append(char)
                    if exceeded:
                        break

                if exceeded:
                    continue

                match.append('*')  # gap symbol

                # After the gap
                for k in range(len_part_2):
                    w_idx = star_index + gap_size + k
                    if w_idx >= len(window):
                        exceeded = True
                        break
                    w = window[w_idx]
                    m = motif_list[star_index + 1 + k]
                    p = penalty_list[star_index + 1 + k]
                    deviation, exceeded, char = check_deviation(w, m, p, deviation, max_deviation)
                    match.append(char)
                    if exceeded:
                        break

                if not exceeded:
                    success = True
                    matched_seq = ''.join(match)
                    break

        if success and deviation <= max_deviation:
            print(f"Yielding match at {i}: {matched_seq}, deviation: {deviation}")
            yield (i, deviation, matched_seq)
        else:
            print(f"No match at {i}, deviation too high: {deviation}")
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
        
        #calculate window length (with or without gap)
        if star_index is not None:
            window = len(motif_list) - 1 + minimum_gap
        else:
            window = len(motif_list)
        
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