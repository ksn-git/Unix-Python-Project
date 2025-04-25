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
                        minimum_gap = int(result.group(1))
                        maximum_gap = int(result.group(2))
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
