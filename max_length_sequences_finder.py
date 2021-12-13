# Given a folder containing fasta formatted files of multiple sequence alignments (the number of sequences per file may vary, but there must be at least one),
# this script finds the file that provides the the longest set of sequences.
# Usage: python3 max_length_sequences.py <folder_path>

import os
import sys

if(len(sys.argv) != 2):
    print("Usage: python3 max_length_sequences.py <folder_path>")
    exit()
else:
    # Set the directory to the folder in question
    os.chdir(sys.argv[1])

    # Open each file in the directory, read the first sequence, and do a length comparison
    longest_file_name = ""
    longest_sequence_length = -1

    for filename in os.listdir():
        # Open the file
        f = open(filename, "r")

        # Read the label of the first sequence
        f.readline()

        sequence_length = 0
        line = f.readline()

        while not ">" in line:
            sequence_length += len(line)
            line = f.readline() # get next line

        # Check if we found a new longest sequence or this is the first file we are checking
        if longest_sequence_length == -1 or sequence_length > longest_sequence_length:
            longest_sequence_length = sequence_length
            longest_file_name = filename

    print(longest_file_name)
        
