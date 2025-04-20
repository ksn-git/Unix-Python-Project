#!/bin/env python3

## Unittest Project

#libraries
import os
import sys
import pytest

#change sys.path to import helper_module found in 'src'
current_path = os.path.dirname(os.path.abspath(__file__))
src_path = os.path.join(current_path,'..','src')
#check if we need to add 'src' to path
if src_path not in sys.path:
    sys.path.insert(0,src_path)

#import functions
from helper_module import add_to_sys_path,get_data_path
add_to_sys_path('src')

### unit test helper_module 
## correctly build sys.path
# using mock to isolate a piece of the code without dependecies
# sometimes uses tmp_path to create an actual temporary directive to test interactions
# makes the test independant of of the running user

#apply fixture with patches
#replaces os with mock object while testing
from unittest.mock import patch
@pytest.fixture
def mock_os_functions():
    with patch('os.path.abspath') as mock_abspath, patch('os.path.isdir') as mock_isdir:
        yield mock_abspath, mock_isdir

## Notes and explanation for fixture and mock tests
#mock_isdir replaces os.path.isdir in test. If True 
# ensures function is good if directory exists.
#moch_abspath replaces os.path.abspath in test.  
# Ensures function is good if example_file is is temp_dir.
#tmp_path creates temporary dir

#test if sys_to_path works and doesn't duplicate path
def test_sys_to_path_function(mock_os_functions,tmp_path):
    mock_abspath, mock_isdir = mock_os_functions
    
    #create tmp_path to provide an actual tmp_dir to interact with
    temp_dir = tmp_path / "example_dir"
    temp_dir.mkdir()
    #mock filepath in temp_dir
    mock_abspath.return_value = str(temp_dir / 'example_file.py')
    mock_isdir.return_value = True

    #set path
    relative_path = 'example_dir'
    expected_path = str(temp_dir)
    
    #add first path to check if work
    add_to_sys_path(relative_path)
    assert expected_path in sys.path

    #add second path to check if it duplicates path
    add_to_sys_path(relative_path)
    sys.path.appent(expected_path) == 1


#if directory does not exist
def test_dir_not_exist(mock_os_functions):
    mock_abspath, mock_isdir = mock_os_functions

    #generate non-existing dir and set path
    mock_abspath.return_value = '/project_root/example_file.py'     #set mocked os.path.abspath 
    mock_isdir.return_value = False                                 #simulates dir does not exist
    relative_path = 'non_existent_dir'

    #test error
    with pytest.raises(FileNotFoundError):
        add_to_sys_path(relative_path)


#edge cases

#unittest if file exist
def test_nonexisting_file():
    with pytest.raises(Exception):
        #get path
        data_path = get_data_path('dn6.fsa')


### unittest reference_motif_TATAAT function
from reference_motif import reference_motif_TATAAT

def test_missing_penalty():
    with pytest.raises(ValueError):
        m1,p1,gap,m2,p2 = reference_motif_TATAAT('missing_penalty_motif.txt')


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