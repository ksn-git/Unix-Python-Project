#!/usr/bin/env python3
# Exercise 7-8
# Peter's Fasta class and get_data_path 
from functions import get_data_path


class Fasta:
    # Alphabets, class variable
    alphabets = {'dna': "ACGT",
                 'dnax': "ACGTNX",
                 'rna': "ACGU",
                 'IUPACdna': "ACGTNXRYMKSWHBVD",
                 'protein': "ACDEFGHIKLMNPQRSTVWY",
                 'proteinx': "ACDEFGHIKLMNPQRSTVWXY",
                 'IUPACprotein': "ABCDEFGHIKLMNPQRSTVWXYZ" }

    # Instatiation, creates the empty list needed
    def __init__(self, filename=None):
        self.headers = list()
        self.sequences = list()
        self.uniqID = set()
        self.pos = 0
        if filename is not None:
            self.load(filename)
    
    def __len__(self):
        return len(self.headers)
    
    def __iter__(self):
        self.pos = 0
        return self
    
    def __next__(self):
        self.pos += 1
        if self.pos > len(self.headers):
            raise StopIteration
        return self.headers[self.pos-1], self.sequences[self.pos-1]
        
<<<<<<< HEAD
    # O(n*m). Scales with the number of sequences and the number of lines in the sequences. 
    def load(self, filename):
=======
    def load(self, filename,subdir='data'):
>>>>>>> 655c9df23d930d8eaae91f84caf3ff79d06b9f9b
        """Reads a fasta file by given filename and creates a list with headers and a list with sequences"""
        self.headers.clear()
        self.sequences.clear()
        self.uniqID.clear()
        # Ensure fasta file can be opened even if in a different folder
        data_path = get_data_path(filename,subdir=subdir)
        # Open the file, not using try-except as I want the expection to propagate to main program.
        infile = open(data_path, 'r')
        for line in infile:
            if line.startswith('>'):
                self.headers.append(line.rstrip())
                self.sequences.append('')
            elif len(self.sequences) == 0:
                continue        # ignore leading non header in file
            else:
                self.sequences[-1] += ''.join(line.split())
        infile.close()
        # Check that identifiers in headers are unique - a requirement in fasta files
        for header in self.headers:
            ID = header.split()[0][1:]
            if ID in self.uniqID:
                raise ValueError(f"IDs in fasta headers are not unique ({ID})")
            else:
                self.uniqID.add(ID)

    def save(self, filename):
        """Expects a filename. Writes a fasta file."""
        # Ready, now write, no error check on purpose
        outfile = open(filename, "w")
        for i in range(len(self.headers)):
            print(self.headers[i], file=outfile)
            for j in range(0, len(self.sequences[i]), 60):
                print(self.sequences[i][j:j+60], file=outfile)
        outfile.close()

    def _verifyrange(self, start, end):
        if start is None:
            return None, None
        if not isinstance(start, int):
            raise ValueError("Argument is a position and must be an integer")
        if 0 > start:
            start = len(self.headers) + start
        if start >= len(self.headers):
            raise IndexError("Accessing outside range of list")
        if end is None:
            return start, None
        if end == "end":
            end = len(self.headers)
        elif not isinstance(end, int):
            raise ValueError('Argument is a position and must be an integer or "end"')
        if 0 > end:
            end = len(self.headers) + end
        if end > len(self.headers):
            raise IndexError("Accessing outside range of list")
        return start, end

    def content(self, start=None, end=None):
        start, end = self._verifyrange(start, end)
        if start is None:
            return self.headers[:], self.sequences[:]
        elif end is None:
            return self.headers[start], self.sequences[start]
        else:
            return self.headers[start:end], self.sequences[start:end]

    def delete(self, start=None, end=None):
        start, end = self._verifyrange(start, end)
        if start is None:
            self.headers.clear()
            self.sequences.clear()
            self.uniqID.clear()
        elif end is None:
            self.uniqID.discard(self.headers[start].split()[0][1:])
            del self.headers[start]
            del self.sequences[start]
        else:
            for i in range(start, end):
                self.uniqID.discard(self.headers[i].split()[0][1:])
            del self.headers[start:end]
            del self.sequences[start:end]

    def insert(self, myheader, mysequence, position=None):
        position, end = self._verifyrange(position, None)
        if isinstance(myheader, str) and isinstance(mysequence, str):
            # Turn into lists for simpler programming later
            myheader = [myheader]
            mysequence = [mysequence]
        if not isinstance(myheader, list) or not isinstance(mysequence, list):
            raise ValueError("Expected two lists or two strings when inserting")
        if len(myheader) != len(mysequence):
            raise ValueError("Header list and sequence list must be of equal length")
        # Prettify and verify
        uniqset = set()
        for i in range(len(myheader)):
            mysequence[i] = ''.join(mysequence[i].split())
            myheader[i] = myheader[i].strip()
            if not myheader[i].startswith(">"):
                myheader[i] = '>' + myheader[i]
            ID = myheader[i].split()[0][1:]
            if ID in uniqset:
                raise ValueError(f"IDs in fasta headers are not unique ({ID})")
            else:
                uniqset.add(ID)
        if len(self.uniqID.intersection(uniqset)) > 0:
            raise ValueError(f"IDs in fasta headers already exists ({ID})")
        # Finally ready to insert
        self.uniqID.update(uniqset)
        if position is None:
            self.headers.extend(myheader)
            self.sequences.extend(mysequence)
        else:
            self.headers[position:position] = myheader
            self.sequences[position:position] = mysequence

    # O(n). Scales with the number of sequences. Regex is assumed to be O(1).
    def verify(self, alphabet, start=None, end=None):
        import re
        if alphabet not in Fasta.alphabets:
            raise KeyError(f"Unknown alphabet: {alphabet}")
        start, end = self._verifyrange(start, end)
        if start is None:
            start, end = 0, len(self.headers)
        elif end is None:
            end = start + 1
        pattern = "[^{}]".format(Fasta.alphabets[alphabet])
        for i in range(start, end):
            # Check both fully upper case and fully lower case alphabets.
            if re.search(pattern, self.sequences[i]) and re.search(pattern.lower(), self.sequences[i]):
                return False
        return True

    def discard(self, alphabet, start=None, end=None):
        import re
        if alphabet not in Fasta.alphabets:
            raise KeyError(f"Unknown alphabet: {alphabet}")
        start, end = self._verifyrange(start, end)
        if start is None:
            start, end = 0, len(self.headers)
        elif end is None:
            end = start + 1
        pattern = "[^{}]".format(Fasta.alphabets[alphabet])
        for i in range(end-1, start-1, -1):
            # Check both fully upper case and fully lower case alphabets.
            if re.search(pattern, self.sequences[i]) and re.search(pattern.lower(), self.sequences[i]):
                del self.headers[i]
                del self.sequences[i]

    def deletethis(self):
        if self.pos == 0 or self.pos > len(self.headers):
            raise IndexError("Not iterating now")
        self.pos -= 1
        self.delete(self.pos)

    def insertthis(self, header, sequence):
        if self.pos == 0 or self.pos > len(self.headers):
            raise IndexError("Not iterating now")
        offset = len(self.headers)
        self.insert(header, sequence, self.pos-1)
        offset = len(self.headers) - offset
        self.pos += offset

    def verifythis(self, alphabet):
        if self.pos == 0 or self.pos > len(self.headers):
            raise IndexError("Not iterating now")
        return self.verify(alphabet, self.pos-1)

    def discardthis(self, alphabet):
        if not self.verify(alphabet, self.pos-1):
            self.deletethis()