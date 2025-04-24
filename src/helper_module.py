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



