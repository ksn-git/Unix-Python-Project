#!/bin/env python3

### add reference_motif function

import sys, re, os
from peter_fasta_class import Fasta
from helper_module import get_data_path

def reference_motif_TATAAT(filename):
    """Finds reference data and saves it as a motif list, penalty list and a gab."""
    #load
    motif_file = get_data_path(filename)

    #open file
    infile = open(motif_file,'r')
    #empty lists 
    motif1,motif2 = [],[]
    penalty1,penalty2 = [],[]
    gap = None
    reading_state = 'motif1'      #starting state

    for line in infile:
        line = line.strip()
        #find empty lines and comments
        if not line or line.startswith('#'):
            #change state depending on comments
            if 'unimportant bases' in line:
                reading_state = 'gap'
            elif '-10' in line:
                reading_state = 'motif2'
            continue

        #detect unimportant bp
        if line.startswith('*'):
            _,value = line.split('\t')
            match = re.match(r'(\d+)-(\d+)',value)
            if match:
                gap = (int(match.group(1)),int(match.group(2)))
            else:
                raise ValueError(f'Invalid gap: {value}')
        
        #detect bp and mismatch score
        content = line.split('\t')
        if len(content) != 2:
            raise ValueError(f'Expected tab-seperated line with 2 columns but got: {line}')
        bp,penalty = content

        #update motif
        if reading_state == 'motif1':
            motif1.append(bp)
            penalty1.append(penalty)
        elif reading_state == 'motif2':
            motif2.append(bp)
            penalty2.append(penalty)
    infile.close()

    return motif1,penalty1,gap,motif2,penalty2


"""
#check
m1,p1,gap,m2,p2 = reference_motif('reference.txt')

print("Motif 1:", m1)
print("Penalty 1:", p1)
print("Gap:", gap)
print("Motif 2:", m2)
print("Penalty 2:", p2)
"""
#unittest needs: what if no gap, no motif after, no motif before,
# no penalty score, no file, totally wrong file type, gap wrong way around

#detect bp and mismatch score found with unittest





