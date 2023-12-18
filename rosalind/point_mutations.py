#!/usr/bin/env python

# Author: Ian VanGordon
# 12/08/2023

'''
This code calculates the hamming distance between two sets of DNA or RNA. The program takes a txt file with each sequence on a different line.
This problem is for "Counting Point Mutations" on Rosalind.
'''

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="This program calculates the hamming distance between two DNA sequences (number of mismatched bps).")
    parser.add_argument("-f", "--file", help="What is the filepath for the txt file to be read?", type=str)
    parser.add_argument("-o", "--outfile", help="What should the output file be called?", type=str)
    return parser.parse_args()

args = get_args()

input_file = args.file
output_file = args.outfile

def read_seq(file):
    seqs = []
    with open(file) as fh:
        for line in fh:
            seqs.append(line.strip())
    return seqs


if __name__ =="__main__":

    seqs = read_seq(input_file)

    seq1 = seqs[0]
    seq2 = seqs[1]
    ham_dist = 0

    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            ham_dist += 1
        else:
            continue

    with open(output_file, "w") as fh:
        fh.write(f"The hamming distance between the two sequences in {input_file} is: {ham_dist}")
