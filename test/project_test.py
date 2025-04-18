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
from peter_fasta_class import Fasta

### unit test helper_module 
## correctly build path
# using mock to isolate a piece of the code without dependecies
# monkey patch fixture safely modifies "object" for importing and 
# makes the test independant of of the running user

#function in isolation

#ensure directory exist

#avoid path duplication 

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
        m1,p1,gap,m2,p2 = reference_motif_TATAAT('missing_penalty_reference.txt')


#unittest needs: what if no gap, no motif after, no motif before,
# no penalty score, no file, totally wrong file type, gap wrong way around


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

