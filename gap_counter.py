# Given a folder containing fasta formatted files of multiple sequence alignments (the number of sequences per file may vary, but there must be at least one),
# this script counts the number of gap in each sequence and prints a list of the labels in descending order of number of gaps, as well as the average number of gaps.                                           
# Usage: python3 gap_counter.py <file_path>

import os
import sys

# count how many time this value occurs in the string
def count(value, string):
    count = 0
    for i in range(0,len(string)):
        if string[i] == value:
            count += 1
    
    return count

# Check usage
if(len(sys.argv) != 2):
    print("Usage: python3 gap_counter.py <file_path>")
    exit()
else:
    filename = sys.argv[1]
        
    # Open the file
    f = open(filename, "r")
    label = ""
    map = {}

    # Read the label of the first sequence
    for line in f:
        if ">" in line: # This is a label
            label = line.strip()[1:]
            map[label] = 0
        else:
            map[label] = map[label] + count("-", line)
    

    result = []
    sum = 0

    while len(result) < len(map):
        max = -1
        max_v = ""

        for key in map:
            if max == -1 or map[key] > max:
                max = map[key]
                max_v = key
        result.append(max_v + ":" + str(map[max_v]))
        sum += map[max_v]
        map[max_v] = -1
    
    print("Sequences in descending order of gap count: " + str(result))
    print("Average Gap Count: " + str(sum/len(map)))
        
