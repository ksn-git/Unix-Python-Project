#!/bin/env python3

## Unittest Project

#libraries
import os
import sys
import pytest
from unittest.mock import patch

#change sys.path to import helper_module found in 'src'
current_path = os.path.dirname(os.path.abspath(__file__))
src_path = os.path.join(current_path,'..','src')
#check if we need to add 'src' to path
if src_path not in sys.path:
    sys.path.insert(0,src_path)

#import functions
from helper_module import add_to_sys_path,get_data_path

### test sys_to_path from helper_module 
# using mock to isolate a piece of the code without dependecies
# sometimes uses tmp_path to create an actual temporary directive to test interactions
# makes the test independant of of the running user

## Notes and explanation for fixture and mock tests
#mock_isdir replaces os.path.isdir in test. If True 
# ensures function is good if directory exists.
#moch_abspath replaces os.path.abspath in test.  
# Ensures function is good if example_file is is temp_dir.
#tmp_path creates temporary dir

#test if sys_to_path works and doesn't duplicate path
def test_sys_to_path_function(tmp_path):
    
    #create tmp_path to provide an actual tmp_dir to interact with
    temp_dir = tmp_path / "example_dir"
    temp_dir.mkdir()

    #mock system from run_from.py
    with patch('os.path.abspath',return_value = str(temp_dir / 'run_from.py')):
        with patch('os.path.isdir',return_value = True):
            #set path
            relative_path = 'example_dir'
            expected_path = str(temp_dir)

            #add first path to check if work
            add_to_sys_path(relative_path)
            assert expected_path in sys.path

            #add second path to check if it duplicates path
            add_to_sys_path(relative_path)
            assert sys.path.count(expected_path) == 1

#if directory does not exist
def test_dir_not_exist():
    with patch('os.path.abspath',return_value = '/project_root'):
        with patch('os.path.isdir',return_value = False):
            relative_path = 'non_existent_dir'
            #tests if FileNotFoundError is detected
            with pytest.raises(FileNotFoundError):
                add_to_sys_path(relative_path)

##edge cases
#check if an empty path is handled correctly
# by raising an exception or ignoring input
def test_sys_path_empty():
    #mock system 
    with patch('os.path.abspath',return_value = '/project_root'):             
        with patch('os.path.isdir',return_value = True):
            #tests if raises exception without relative_path
            with pytest.raises(ValueError):
                add_to_sys_path('')

### test get_data_path from helper_module
#test correct path
def test_get_correct_data_path(tmp_path):
    #create temp file location
    filename = 'testfile.txt'
    temp_dir = tmp_path / 'data'
    temp_dir.mkdir()
    expected_path = temp_dir / filename
    #create file
    expected_path.write_text('content')

    #mock system from run_from.py
    with patch('os.path.abspath',return_value = str(temp_dir / 'run_from.py')):
        assert get_data_path(filename) == str(expected_path)

#test file does not exist
def test_get_data_path_file_not_found(tmp_path):
    #don't create file
    filename = 'non_existing_file.txt'

    #mock system from run_from.py
    with patch('os.path.abspath',return_value = str(tmp_path / 'run_from.py')):
        with pytest.raises(FileNotFoundError) as err:
            get_data_path(filename)



#unittest needs: what if no gap, no motif after, no motif before,
# no penalty score, no file, totally wrong file type, gap wrong way around

"""
### find motif generator
from main_program import find_motif

#Basic Functionality Test
def test_find_motif_basic():
    sequence = 'CGCCTATAATAAT'
    motif = ['T','A','T','A','A','T']
    penalty = [8,8,6,6,5,8]
    max_deviation = 0
    result = list(find_motif(sequence,motif,penalty,max_deviation))
    assert result == [(4,0,'TATAAT')]

#functionality test with gap
def test_find_motif_with_gap():
    sequence = 'ATCGGACCCACTAGT'
    motif = ['A','T','C','G','G','A','*','*','*','A','G','T','C','G','T']
    penalty = [8,7,6,7,8,9,0,0,0,8,7,6,9,8,4] 
    max_deviation = 18
    result = list(find_motif(sequence,motif,penalty,max_deviation))
    assert result == [(0,16,'ATCGGACCCACTAGT')]

#with multiple available bp
def test_find_motif_multiple_bp():
    sequence = 'CGCCTATAATAAT'
    motif = ['T','A','T','A','AT','T']
    penalty = [8,8,6,6,5,8]
    max_deviation = 0
    result = list(find_motif(sequence,motif,penalty,max_deviation))
    assert result == [(4,0,'TATAAT')]

#Penalty Application Test (need testing)
def test_find_motif_penalty():
    sequence = 'CGCCTATCATTATCCT'
    motif = ['T','A','T','A','A','T']
    penalty = [8,8,6,6,5,8]
    max_deviation = 8
    result = list(find_motif(sequence,motif,penalty,max_deviation))
    assert result == [(4,8,'TATCAT')]

## Edge Case Test
#empty sequence
def test_find_motif_empty_sequence():
    sequence = ''
    motif = ['T','A','T','A','A','T']
    penalty = [8,8,6,6,5,8]
    max_deviation = 0
    result = list(find_motif(sequence,motif,penalty,max_deviation))
    assert result == []

#empty motif and penalty
def test_find_motif_empty_motif_and_penalty():
    sequence = 'CGCCTATAATAAT'
    motif = []
    penalty = []
    max_deviation = 0
    result = list(find_motif(sequence,motif,penalty,max_deviation))
    assert result == []

#empty penalty score
def test_find_motif_empty_motif_and_penalty():
    sequence = 'CGCCTATAATAAT'
    motif = ['T','A','T','A','A','T']
    penalty = []
    max_deviation = 0
    result = list(find_motif(sequence,motif,penalty,max_deviation))
    assert result == []

#ensure fasta file
#meet another star
"""
