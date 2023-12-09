#!/usr/bin/env python

# Author: Ian VanGordon
# 12/08/2023

'''
This code finds all reverse palindromes or restriction sites from a sequence of DNA.
This script solves the Rosalind challenge "Locating Restriction Sites"
'''

import argparse
import sys
sys.path.append("..")
from bioinfo_toolbox import DNA_compliments

def get_args():
    parser = argparse.ArgumentParser(description="This program locates reverse palindromes in a DNA sequence.")
    parser.add_argument("-f", "--file", help="What is the filepath for the txt file to be read?", type=str)
    parser.add_argument("-o", "--outfile", help="What should this file be called?", type=str)
    return parser.parse_args()


def find_rev_pal(fw_seq):
    '''This function takes a sequence and returns the position and length of every reverse palindrome between 4 and 12 bp long.'''
    pos_length = []
    for i in range(len(fw_seq)):                                                                                # Iterates through each position in the sequence
        for j in range(i, (len(fw_seq))):                                                                       # Compares every bp foreward to the current bp in the loop above
            dist = j - i + 1                                                                                    # Distance between the bps in the two loops
            if (fw_seq[i] == DNA_compliments[fw_seq[j]]) and (dist > 2) and (dist % 2 == 0) and (dist <= 12):   # Finds complimentary nucleotides and ensure that the distance is even and between 4 and 12 bp away
                half_dist = int(dist / 2)
                counter = 0
                for k in range(half_dist):                                                                      # Iterates the two current base pairs getting closer to the middle making sure that they are reverse palindromes
                    if fw_seq[i + k] == DNA_compliments[fw_seq[j - k]]:
                        counter += 1
                if counter == half_dist:                                                                        # Adds position and distance to the list of tuples to be called on later
                    cur = ((i + 1), dist)
                    pos_length.append(cur)
    return pos_length

if __name__ == "__main__":
    args = get_args()

    input_file = args.file
    output_file = args.outfile

    counter = 0
    fw_seq = ""

    with open(input_file) as fh:
        for line in fh:
            if line[0] != ">":
                fw_seq += line.strip()

    pos_length = find_rev_pal(fw_seq)

    with open(output_file, "w") as fh:
        for i in range(len(pos_length)):
            fh.write(f"{pos_length[i][0]}\t{pos_length[i][1]}\n")

