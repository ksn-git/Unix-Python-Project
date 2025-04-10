#!/bin/env python3

## Unittest Project

# Secure path
import os
import sys
#get parent directive (of test and data)
project_path = os.path.abspath(os.path.join(os.path.dirname(__file__),'..')) 
#add testdata to directive  
src_path = os.path.join(project_path,'src')                 
sys.path.insert(0,src_path)
#import from new path 
from peter_fasta_class import Fasta

#module for unit test
import pytest

### load 
#check correct responses
@pytest.mark.parametrize('filename',['dna7.fsa','ex1.dat','motif.fsa'])
def test_Fastaloadfile(filename):
    #redefine path to files
    import os
    project_path = os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))   #get parent directive (of test and testdata)
    data_path = os.path.join(project_path,'data',filename)                 #add testdata to directive
    #if file/path exist
    assert os.path.exists(data_path), f'File {filename} does not exist or is placed elsewhere.'

    #run class
    myfasta = Fasta()
    myfasta.load(data_path)

#check IOError response
@pytest.mark.parametrize(('filename'),["ABC"])
def test_fastaload_error(filename):
    with pytest.raises(IOError):
        myfasta = Fasta()
        myfasta.load(filename)

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
    assert "Hello\n" == txt

"""
def save(self, filename):
        """Expects a filename. Writes a fasta file."""
        # Ready, now write, no error check on purpose
        outfile = open(filename, "w")
        for i in range(len(self.headers)):
            print(self.headers[i], file=outfile)
            for j in range(0, len(self.sequences[i]), 60):
                print(self.sequences[i][j:j+60], file=outfile)
        outfile.close()"
"""
#have headers
#same number of headers and seq
#does content fit?


