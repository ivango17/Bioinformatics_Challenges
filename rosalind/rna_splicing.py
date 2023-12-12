#!/usr/bin/env python

# Author: Ian VanGordon
# 12/11/2023

'''
This script takes a FASTA file with the first sequence being a whole region of DNA and all of the subsequent sequences are introns to be spliced out. The the resulting protein string is output.
The program solves the Rosalind challenge RNA Splicing.
'''

import argparse
import re
import sys
sys.path.append("..")
from bioinfo_toolbox import dna_to_aa


def get_args():
    parser = argparse.ArgumentParser(description="This program takes a DNA sequence with introns and splices the DNA to then return the protein coded in that sequence.")
    parser.add_argument("-f", "--file", help="What is the filepath for the txt file to be read?", type=str)
    return parser.parse_args()


def rna_splice(seq, intron_list):
    '''This function takes a string of RNA and a list of introns. It returns a string of spliced RNA.'''
    for intron in intron_list:
        if re.split(intron, seq):
            exon_list = re.split(intron, seq)
            seq = "".join(exon_list)
            
    return seq


def prot_string(seq):
    '''This function checks a protein string for stop codons and removes the stop codon at the end of the string. If there are other stop codons in the string the program will exit.'''
    if "*" in seq : seq = seq[0:-1]
    if "*" in seq : print("Protein string contains multiple stop codons!"); exit()
    return seq


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

    seq = fasta_seqs.pop(0)

    result = rna_splice(seq, fasta_seqs)
    result = dna_to_aa(result)
    result = prot_string(result)
    print(result)





