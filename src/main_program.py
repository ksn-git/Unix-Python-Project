#!/bin/env python3

### Get imports and arguments from commandline
#example of use with our files on cmd
# ./main_program.py motif.txt sequence.fsa 20 data matched_seq.txt
# python3 main_program.py motif.txt sequence.fsa 20 data matched_seq.txt

# Dependencies
import sys
from peter_fasta_class import Fasta
from functions import get_data_path,load_motif,find_motif

# Check length of arguments
if len(sys.argv) != 6:
    print("You must supply all needed inputs. Either:")
    print("Usage: ./main_program.py <motif> <fastafile> <max deviation> <filedir> <save_as_file>")
    print("or")
    print("Usage: python3 main_program.py <motif> <fastafile> <max deviation> <filedir> <save_as_file>")
    sys.exit(1)

# Unpack arguments 
motif,fasta_file,max_deviation,subdir,save_as_file = sys.argv[1:]                              

# Check and update max_deviation data type
try:
    max_deviation = int(max_deviation)
except ValueError:
    print(f'max_deviation must an integer and got: {max_deviation}')
    sys.exit(1)

### Usage of code 
# Loading and verifying the fasta file 
fasta = Fasta()
fasta.load(fasta_file,subdir=subdir)   
fasta.verify("dnax")     

# Opening the motif description file. Works if the file is in subdir argument or default:"data" subfolder. 
motif_file = get_data_path(motif,subdir=subdir)       
# Load motif and save the function output
motif_list, penalty_list, minimum_gap, maximum_gap = load_motif(motif_file)
#print(motif_list, penalty_list, minimum_gap, maximum_gap)      #debug statement

#add matches into dict
count = 0       #debug use
output = {}
for header,sequence in fasta:
    #add new fasta headers
    if header not in output:
        output[header] = []
    #append matches
    for match in find_motif(sequence, motif_list, penalty_list, max_deviation, minimum_gap, maximum_gap):
        output[header].append(match)
        count += 1  #debug use

#read headers and matches into file
outputfile = get_data_path(save_as_file,subdir,must_exist = False)
with open(outputfile,'w') as outfile:
    #add comment(s) to file
    outfile.write("# Matches are listed as (start position, penalty score, match)\n")
    outfile.write("# The header corresponding to the match is printed immediately before the match\n")
    #print to terminal
    print("Matches are listed as (start position, penalty score, match)")
    print("The header corresponding to the match is printed immediately before the match")
    #add content
    for header,matches in output.items():
        outfile.write(header + '\n')    #add to outfile
        print(header)                   #print to terminal
        for pos,dev,mstr in matches:
            outfile.write(f'{pos}\t{dev}\t{mstr}\n')    #add to outfile
            print(f'({pos},{dev},{mstr})')              #print to terminal
print(f'Output has been saved in: {save_as_file}')

#debug statement to see number of matches
print(count)
