#!/usr/bin/env python

# Author: Ian VanGordon
# 12/11/2023

'''
This script takes a FASTA file with two sequences and calculates the transition to transversion ratio between the two sequences.
This program solves the Rosalind challenge Transitions and Transversions.
'''

import argparse
import sys
sys.path.append("..")
from bioinfo_toolbox import transversions, transition_transverion

def get_args():
    parser = argparse.ArgumentParser(description="This program takes a FASTA file of DNA sequences and returns the longest motif or motifs that are shared between all of the reads.")
    parser.add_argument("-f", "--file", help="What is the filepath for the txt file to be read?", type=str)
    return parser.parse_args()


if __name__ == "__main__":

    args = get_args()

    input_file = args.file

    cur_seq = ""
    fasta_seqs = []

    with open(input_file) as fh:
        for line in fh:
            if line[0] != ">":
                cur_seq += line.strip()
            elif cur_seq != "":
                fasta_seqs.append(cur_seq)
                cur_seq = ""
        fasta_seqs.append(cur_seq)

    print(transition_transverion(fasta_seqs[0], fasta_seqs[1]))