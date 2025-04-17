
import sys, re, os
from peter_fasta_class import Fasta
from helper_module import get_data_path

# Getting filename with motif from command line
if len(sys.argv) != 4:                              
    if len(sys.argv) == 1:
        motif = input("Please enter the name of the file with the motif to look for: ")
        fasta_file = input("Please enter the name of the fasta file to search in: ")
        max_deviation = int(input("Please enter the maximum deviation (int): "))
    else:
        print("You must supply a file name with a motif to look for")
        print("Usage: <main_program.py> <motif> <fastafile> <max deviation>")
        sys.exit(1)
else:
    motif = sys.argv[1]                              
    fasta_file = sys.argv[2]
    max_deviation = int(sys.argv[3])

# Loading and verifying the fasta file 
fasta = Fasta()
fasta.load(fasta_file)                           
fasta.verify("dnax")     

# Opening the motif description file. Works if the file is in a "data" subfolder. 
# check if file exist is in helper_module 
motif_file = get_data_path('motif.fsa')     #this is the fasta file, I've added a reference file     
infile = open(motif_file,'r')

# Arranging motif into two lists
# Motif list:   ATCGGATC******AGTCGTTA
# Penalty list: 8767894200000087698432
motif_list = []
penalty_list = []
for line in infile:
    # Skip descriptive lines
    if line.startswith("#"):
        continue
    else:
        # Line specifying unimportant positions
        if line.startswith("*"):
            unimportant_positions = line.strip().split(sep="\t")[1]
            # If there is a defined number of unimportant positions
            try:
                unimportant_positions = int(unimportant_positions)
                for i in range(0, unimportant_positions):
                    motif_list.append("*")
                    penalty_list.append(0)
                    continue 
            except ValueError:
                pass                
            # If there is a range of unimportant positions
            if isinstance(unimportant_positions, str):
                result = re.search(r"(\d+)-(\d+)", unimportant_positions)
                if result is not None:
                    minimum = int(result.group(1))
                    maximum = int(result.group(2))
                    for i in range(0, minimum):
                        motif_list.append("*")
                        penalty_list.append(0)
                else:
                    raise ValueError("Invalid unimportant positions format")
                    sys.exit(1)
        # Lines specifying important positions
        else:
            motif_list.append(line.strip().split(sep="\t")[0])
            penalty_list.append(line.strip().split(sep="\t")[1])
infile.close()

#find motif
#I think this should be added to the class...
def find_motif(sequence,motif_list,penalty_list,max_deviation):
    """Generator that yields a matching motif."""
    #to avoid index error
    print('placeholder for now')

    ##psudocode
    
    # for seq in sequences
        #reset mismatch score
        #reset temp match 

        # for each sequence position in seq (using sliding window)
            #check if sequence position value is below deviation
            # if temp match != '' and last postion in match is reached (x positions filled with not '*'?)
                #save/yield? (temp match,position,mismatch score) to match list for the first entry
                #reset mismatch score
                #reset temp match
            # if over deviation,
                #move to next start position
                #reset mismatch score
                #reset temp match
            # if below deviation,
                #add possible mismatch value to mismatch score
                #add position to temp match
            #if match is in progress and skip position is reached (start to 15)
                #add '*' to file and don't 
    
# need to find out how to jump back to last starting position. 
# might be missing some conditions


## commented out for now 
'''
def find_motif(sequence, motif_list, penalty_list, max_deviation):
    """ Generator that yields a result when a motif is found"""

    # Make a set for each reading frame. Elements of 3 bases from motif. 
    reading_frame1 = list()
    reading_frame2 = list()
    reading_frame3 = list()

    # OBS: ONLY WORKS WHEN EACH POSITION HAS ONE POSSIBLE LETTER
    for i in range(0, len(motif_list) - 2):
        codon = str(motif_list[i] + str(motif_list[i + 1]) + str(motif_list[i + 2]))

        # Include codons with *. This will make tracking codon position easier later.
        # Reading frame 1
        if i % 3 == 0:
            reading_frame1.append(codon)
        # Reading frame 2
        if (i + 2) % 3 == 0:
            reading_frame2.append(codon)
        # Reading frame 3
        if (i + 1) % 3 == 0:
            reading_frame3.append(codon)

    # Sliding window of 3 bases through the sequence. Looking for codons in the reading frames.
    pos_in_rf1_list = list()
    pos_in_rf2_list = list() 
    pos_in_rf3_list = list()      
    start_pos_rf1 = list()
    start_pos_rf2 = list()
    start_pos_rf3 = list()
    for i in range(0, len(sequence) - 2):
        window = sequence[i:i + 3]
        # Reading frame 1
        if i % 3 == 0:
            if window in reading_frame1:
                pos_in_rf1_list.append(reading_frame1.index(window))     # Which # codon in the reading frame was it? 
                start_pos_rf1.append(i)                                  # Where in the sequence was it found?
        # Reading frame 2
        if (i + 2) % 3 == 0:
            if window in reading_frame2:
                pos_in_rf2_list.append(reading_frame2.index(window))     
                start_pos_rf2.append(i)
        # Reading frame 3
        if (i + 1) % 3 == 0:
            if window in reading_frame3:
                pos_in_rf3_list.append(reading_frame3.index(window))     
                start_pos_rf3.append(i)

find_motif(fasta.sequences[0], motif_list, penalty_list, max_deviation)
'''

# Searching for the motif in each entry in the fasta file
#for sequence in fasta.sequences:
#    find_motif(motif_list, penalty_list, sequence, max_deviation)


