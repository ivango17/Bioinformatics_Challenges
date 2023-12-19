#!/usr/bin/env python

# Author: Ian VanGordon
# 12/18/2023

'''
This program takes a FASTA sequence and returns the kmer composition as frequencies in lexicongraphic order.
This script solves the Rosalind challenge k-Mer Composition.
'''

import argparse
from enumerating_kmers import combinations
import sys
sys.path.append("..")
from bioinfo_toolbox import kmerize

def get_args():
    parser = argparse.ArgumentParser(description="This program takes a string of space separated characters and returns every combination of r length of those characters in lexicographic order.")
    parser.add_argument("-f", "--file", help="Where is the file containing the FASTA sequence?", type=str)
    parser.add_argument("-o", "--outputfile", help="What is the destination for the output?", type=str)
    return parser.parse_args()


def kmer_frequency(kmer_list, poss_kmers):
    '''This function takes a list of kmers length k and a list of all possible kmers of length k. It returns the kmer composion frequencies in lexicographic order.'''
    kmer_freq = {}
    lex_freq_order = []
    for kmer in kmer_list:
        if kmer in kmer_freq:
            kmer_freq[kmer] += 1
        else:
            kmer_freq[kmer] = 1

    for poss in poss_kmers:
        if poss in kmer_freq:
            lex_freq_order.append(kmer_freq[poss])
        else:
            lex_freq_order.append(0)

    return lex_freq_order


if __name__ == "__main__":

    args = get_args()

    file = args.file
    output = args.outputfile

    bases = ['A', 'C', 'G', 'T']
    seq = ""

    with open(file) as fh:
        for line in fh:
            if line[0] != ">":
                seq += line.strip()


    kmer_list = sorted(kmerize(seq, 4))
    poss_kmers = combinations(bases, 4)

    kmer_comp = kmer_frequency(kmer_list, poss_kmers)

    with open(output, "w") as fh:
        for freq in kmer_comp:
            fh.write(f"{freq} ")

