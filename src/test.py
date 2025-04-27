#!/bin/env python3

#dependency load_motif
import re

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
                minimum_gap, maximum_gap = gap_positions, gap_positions
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
            penalty_list.append(row[1])
        # Lines specifying important positions   
        else:
            motif_list.append(row[0])
            penalty_list.append(row[1])
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
        len_part_2 = 0      

    #find length of input
    seq_len = len(sequence)
    motif_len = len(motif_list) - (1 if star_index is not None else 0)          #substract 1 if gap exist
    max_window = motif_len + (maximum_gap if star_index is not None else 0)     #add gap to window if it exists

    # Start searching for matches
    for i in range(len(seq_len - max_window + 1)):                            # Search until the remaining seq is not long enough to be the motif
        # reset match
        window = sequence[i:i+max_window]
        deviation = 0
        index = 0                               #index in window or position?
        while index < len(window):
            # while the gap exist and the gap is reached
            if star_index is not None and index == star_index:
                #try with all possible gaps
                for gap_position in range(minimum_gap,maximum_gap +1):
                    deviation_gap = 0      
                    valid = True
                    for m in range(len_part_2):
                        window_pos = index + gap_position + m 
                        motif_pos = star_index + 1 + m 
                        #break if exceeding window or motif
                        if window_pos >= len(window) or motif_pos >= len(motif_list):
                            valid = False
                            break
                        expected = motif_list[motif_pos]            #match
                        actual = window[window_pos]                 #actual character
                        penalty = penalty_list[motif_pos]           #penalty
                        #if multiple allowed characters are found
                        if isinstance(expected,set):
                            if actual not in expected:
                                deviation_gap += penalty
                        else:
                            if actual != expected:
                                deviation_gap += penalty
                        if deviation_gap > max_deviation:
                            valid = False
                            break
                    if valid and deviation + deviation_gap <= max_deviation:
                        start = star_index + gap_position
                        end = star_index + gap_position + len_part_2
                        matched_sequence = window[:star_index] + window[start:end]
                        yield (i,deviation + deviation_gap,matched_sequence)

                        