#!/bin/env python3

#libraries
import os
import sys

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
    #dependencies
    import re

    infile = open(motif,'r')
    motif_list = []
    penalty_list = []
    for line in infile:
        # Skip descriptive lines
        if line.startswith("#"):
            continue
        # Line specifying unimportant positions
        elif line.startswith("*"):
            unimportant_positions = line.strip().split(sep="\t")[1]
            # If there is a defined number of unimportant positions
            try:
                unimportant_positions = int(unimportant_positions)
                minimum = unimportant_positions
                maximum = unimportant_positions
                motif_list.append("*")
                penalty_list.append(0) 
            except ValueError:
                pass                
            # If there is a range of unimportant positions
            if isinstance(unimportant_positions, str):
                result = re.search(r"(\d+)-(\d+)", unimportant_positions)
                if result is not None:
                    minimum = int(result.group(1))
                    maximum = int(result.group(2))
                    motif_list.append("*")
                    penalty_list.append(0)
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