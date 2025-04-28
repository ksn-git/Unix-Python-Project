#!/bin/env python3

#check deviation
def check_deviation(window_elem,motif_elem,penalty_elem,deviation,max_deviation):
    """Compare the element from the window with the motif element to return the updated deviation
    and boolean indicating success or not (exceeding max)"""

    #if multiple charactors are found
    if isinstance(motif_elem,set):
        if window_elem not in motif_elem:
            deviation += penalty_elem
            if deviation > max_deviation:
                return deviation, True      #signal stop now
    #if a single charactor is found
    else: 
        if window_elem != motif_elem:
            deviation += penalty_elem
            if deviation > max_deviation:
                return deviation, True
    #return False if match continues
    return deviation, False         #continue match


#find motif - not fully functional yet!
def find_motif(sequence, motif_list, penalty_list, max_deviation, minimum_gap, maximum_gap):
    """ Generator that searches for motif """
    """ Yields position, deviation and sequence when a match is found """
    """ Tries to make a match in each window until max deviation is reached """
    # Find start position of gap and length of 2nd part of the motif
    try:
        star_index = motif_list.index("*")      
        len_part_2 = len(motif_list) - star_index - 1
        print(f"Gap (star) found at index {star_index}, length of second part of motif: {len_part_2}")
    # If no gaps in motif, ignore
    except ValueError:
        star_index = None
        print("No gap (star) found in motif.")

    # Start searching for matches
    for i in range(0, len(sequence) - len(motif_list) - maximum_gap + 1):           # Search until the remaining seq is not long enough to be the motif
        print(f"\nChecking window starting at position {i}")
        # If no gaps in the motif, window is the same length as the motif
        if minimum_gap == 0 and maximum_gap == 0:
            window = sequence[i:(i + len(motif_list))]                
        # If there are gaps in the motif, window is the motif length + maximum number of gaps in motif. 
        # Subtract 1 because there is a star character in the motif list denoting gaps
        else:
            window = sequence[i:(i + len(motif_list) - 1 + maximum_gap)]                # Window of the same length as the longest possible motif. Len = 29
        print(f"Window: {window}")

        # Reset deviation score and search window for motif
        deviation = 0
        for j in range(len(window)):
            print(f"Checking position {j}, motif element: {motif_list[j]} vs window element: {window[j]}")

            # If gap has been reached 
            # Check match for all gap sizes. Yield if below max deviation
            if j == star_index:
                # Range of gaps
                for k in range(minimum_gap, maximum_gap + 1): # range should be the length of the 2nd part of the motif
                    print(f"Gap at position {j}, checking for gap size {k}")

                    # check match for the length of the 2nd part of the motif
                    for m in range(k, k+len_part_2): # maybe minus 1
                        print(f"Checking position {m} in window for match: {window[j+m]} vs motif element: {motif_list[j+m-k+1]}")

                        # check match and discontinue if exceeding max_deviation
                        deviation, max_exceeded = check_deviation(window[j+m],motif_list[j+m-k+1], penalty_list[j+m-k+1], deviation, max_deviation)
                        print("Max deviation exceeded, breaking out of gap checking loop.")

                        if max_exceeded:
                            break  
                    
            # If gap has been reached, the entire window has already been checked
            elif star_index is not None and j >= star_index:
                # Skip this part if the star/gap has already been encountered
                continue
            
            else:
                # check match and discontinue if exceeding max_deviation
                deviation, max_exceeded = check_deviation(window[j],motif_list[j], penalty_list[j], deviation, max_deviation)
                if max_exceeded:
                    print(f"Max deviation exceeded at position {j}, breaking out of the loop.")

                    break

        # If the deviation is less than the max, yield the match within the window
        if deviation <= max_deviation:
            matched_sequence = window[:star_index] + window[star_index + k:star_index + k + len_part_2]
            yield((i, deviation, matched_sequence))  
