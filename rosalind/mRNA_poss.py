#!/usr/bin/env python

# Author: Ian VanGordon
# 12/10/2023

'''
This script returns the number of possible mRNA sequences that a protein sequence could have been translated from.
This script solves the Rosalind challenge "Inferring mRNA from a Protein."
'''

import argparse
import sys
sys.path.append("..")
from bioinfo_toolbox import aa_numcodon

def get_args():
    parser = argparse.ArgumentParser(description="This program calculates the number of possible mRNA sequences that a protein sequence could have been derived from.")
    parser.add_argument("-f", "--file", help="What is the filepath for the txt file to be read?", type=str)
    return parser.parse_args()

def num_possible_seqs(seq, mod):
    '''This function returns the total number of possible mRNA strands mod some number from a given protein string. '''
    possible = 1
    for base in seq:
        possible *= aa_numcodon[base] 

    possible *= aa_numcodon["*"]
    possible = possible % mod
    return possible


if __name__ == "__main__":

    args = get_args()

    input_file = args.file

    counter = 0
    aa_seq = ""

    with open(input_file) as fh:
        for line in fh:
            if line[0] != ">":
                aa_seq += line.strip()

    print(num_possible_seqs(aa_seq, 1_000_000))


