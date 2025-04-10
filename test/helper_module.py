#!/bin/env python3

#libraries
import os
import sys

#change sys_path
def add_to_sys_path(relative_path_from_root: str):
    """Add directory to sys path from project root."""

    #find current path
    current_file = os.path.abspath(__file__)
    #go back and look at project root
    project_root = os.path.dirname(os.path.dirname(current_file))
    #get target path
    target_path = os.path.join(project_root,relative_path_from_root)

    #add target_path if needed
    if target_path not in sys.path:
        sys.path.insert(0,target_path)
    
