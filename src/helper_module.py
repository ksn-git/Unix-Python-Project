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


# detect motif
def detect_motif(motif):
    """Loads and detects motif from file.
    It arranges the motif into two lists (example):
     Motif list:   ['A','T','C','G','G','A','*','*','*','A','G','T','C','G','T']
     Penalty list: [8,7,6,7,8,9,0,0,0,8,7,6,9,8,4] 
     """

    infile = open(motif,'r')
    motif_list = []
    penalty_list = []
    for line in infile:
        # Skip descriptive lines
        if line.startswith("#"):
            continue
        #seperate line into parts
        parts = line.strip().split('\t')
        symbol = parts[0]

        # Line specifying unimportant positions
        if symbol == '*':
            num_positions = parts[1]
            # If a fixed number of unimportant positions are found (ex. 15)
            try:
                count = int(num_positions)      
                motif_list.extend(['*'] * count)
                penalty_list.append([0] * count) 
            except ValueError:              
                # check if a range of unimportant positions are found (ex. 15-21)
                match = re.search(r"(\d+)-(\d+)", num_positions)
                if match:
                    minimum = int(match.group(1))
                    maximum = int(match.group(2))
                    count = maximum                 #assume maximum number of spaces
                    motif_list.append(['*'] * count)
                    penalty_list.append([0] * count)
                else:
                    raise ValueError("Invalid unimportant positions format")
        
        # Lines with >1 possible character 
        elif len(symbol) > 1:
            motif_list.append(set(symbol))
            penalty_list.append(int(parts[1]))
        # Lines specifying important positions   
        else:
            motif_list.append(symbol)
            penalty_list.append(int(parts[1]))
    infile.close()

    return motif_list,penalty_list
