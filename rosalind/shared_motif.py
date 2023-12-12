#!/usr/bin/env python

# Author: Ian VanGordon
# 12/10/2023

'''
This program takes a FASTA file of DNA sequences and returns the longest motif or motifs that are shared between all of the reads.
This program solves the Rosalind challenge Finding a Shared Motif.
'''

import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="This program takes a FASTA file of DNA sequences and returns the longest motif or motifs that are shared between all of the reads.")
    parser.add_argument("-f", "--file", help="What is the filepath for the txt file to be read?", type=str)
    return parser.parse_args()


def find_motif(seq1, seq2, max_len):
    '''This function compares two sequences and returns the longest motif(s).'''
    print(len(seq1), " ", len(seq2))
    if len(seq1) != len(seq2):
        lg_seq = max(seq1, seq2)
        sm_seq = min(seq1, seq2)
    else:
        lg_seq = seq1
        sm_seq = seq2
    i = min(len(sm_seq), max_len)
    if sm_seq not in lg_seq:
        while i >= 2:
            i -= 1
            kmer_list = [sm_seq[kmer:kmer + i] for kmer in range(0, len(sm_seq) - (i-1))]
            shared_motif = [motif for motif in kmer_list if re.search(motif, lg_seq)]
            
            if shared_motif:
                return shared_motif
    else:
        return [sm_seq]


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

    motif = [fasta_seqs[0]]
    max_motif = []

    for i in range(1, len(fasta_seqs)):                             # This loop could be cleaned up

        for j in range(len(motif)):
            temp_motif = find_motif(motif[j], fasta_seqs[i], 300)
            if temp_motif == None:
                continue
            elif j == 0 and (temp_motif):
                max_motif = temp_motif[0:len(temp_motif)]

            else:
                if (max_motif) and (temp_motif):
                    print("executing")
                    if max(temp_motif) > max(max_motif):
                        max_motif = temp_motif[0:len(temp_motif)]
                else:
                    max_motif = temp_motif[0:len(temp_motif)]
        motif = max_motif[0:len(max_motif)]



    print(motif)