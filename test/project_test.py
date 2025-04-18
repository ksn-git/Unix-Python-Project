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

# it adds path if not already present

# it doesn't add the path multiple times

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








### load 
#check correct responses has sequences with content
@pytest.mark.parametrize('filename',['dna7.fsa','motif.fsa'])
def test_Fastaloadfile(filename):
    #get full path
    data_path = get_data_path(filename)
    #load Fasta class
    myfasta = Fasta()
    myfasta.load(data_path)
    
    #check load
    assert hasattr(myfasta,'sequences'), f'No sequences attribute exist in {filename}.'
    assert len(myfasta.sequences) > 0, f'Sequnece list is empty in {filename}'


#check IOError response
@pytest.mark.parametrize(('filename'),["ABC"])
def test_fastaload_error(filename):
    with pytest.raises(IOError):
        myfasta = Fasta()
        myfasta.load(filename)


"""
#check if the content of the loaded file is right
@pytest.mark.parametrize('filename',['dna7.fsa','motif.fsa'])
def test_Fastaload_content(filename):
    #redefine path to files
    import os
    project_path = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))   #get parent directive (of test and testdata)
    data_path = os.path.join(project_path,'data',filename)                 #add testdata to directive
    #if file/path exist
    assert os.path.exists(data_path), f'File {filename} does not exist or is placed elsewhere.'

    #run class
    myfasta = Fasta()
    myfasta.load(data_path)

    #check if headers and sequences exist
    assert len(myfasta.headers) > 0, 'FASTA files needs headers'
    #check if sequences are dna sequences
    dna = {'A','T','C','G','N'}                 #include unknown base N
    for i,seq in enumerate(myfasta.sequences):
        assert set(seq.upper()).issubset(dna), 'Characters outside of DNA was found in the sequences'
    #length comparison (doesn't do anything for now)
    #assert len(myfasta.headers) == len(myfasta.sequences), 'There must be an equal amout of headers and sequences.'

#load needs a more strict check of the sequences (class or when using the sequence)
# test checks for it -> add to function


### save 
#save temp file
def test_writefile(tmp_path):
    filename = tmp_path / "myfile"
    outfile = open(filename, "w")
    print('Hello', file=outfile)
    outfile.close()
    infile = open(filename, "r")
    txt = infile.readline()
    infile.close()
    assert 'Hello\n' == txt

"""

### find motif generator
from main_program import find_motif

#Basic Functionality Test
def test_basic_sequence():
    sequence = 'CGCCTATAATAAT'
    motif = 'TATAAT'
    penalty = [8,8,6,6,5,8]
    max_deviation = 0
    result = list(find_motif(sequence,motif,penalty,max_deviation))
    assert result == [(4,0,'TATAAT')]

#Penalty Application Test (need testing)
def test_penalty():
    sequence = 'CGCCTATACT'
    motif = 'TATAAT'
    penalty = [8,8,6,6,5,8]
    max_deviation = 8
    result = list(find_motif(sequence,motif,penalty,max_deviation))
    assert result == [(4,8,'TATACT')]

#Max Deviation Test


#Edge Case Test

