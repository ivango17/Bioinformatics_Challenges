#!/usr/bin/env python

# Author: Ian VanGordon
# 12/12/2023

'''
This script takes a FASTA file of DNA sequences and returns the adjacency list of of directed edges (tail, head) for matching k overlap.
This program solves the Rosalind challenge Overlap Graphs.
'''

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="This script takes a FASTA file of DNA sequences and returns the adjacency list of of directed edges (tail, head) for matching k overlap.")
    parser.add_argument("-f", "--file", help="What is the filepath for the txt file to be read?", type=str)
    parser.add_argument("-o", "--outfile", help="What should the output file be called?", type=str)
    return parser.parse_args()

def th_overlap(tail, head, overlap):
    '''This function takes two strings and determines if the last three characters of tail match the first three of head.'''
    if tail[len(tail)-overlap:] == head[0:overlap]:
        return True

    else:
        return False


def ol_graph(fasta_dict, overlap):
    '''This function takes a dictionary of reads from a FASTA with the keys being headers and the values being the sequence. It then returns a list of directed edges from the dictionary.'''
    edges = []
    for key1 in fasta_dict:
        for key2 in fasta_dict:
            if (key1 != key2) and (fasta_dict[key1] != fasta_dict[key2]):
                if (th_overlap(fasta_dict[key1], fasta_dict[key2], overlap)) and ((key1, key2) not in edges):
                    edges.append((key1, key2))
                if (th_overlap(fasta_dict[key2], fasta_dict[key1], overlap)) and ((key2, key1) not in edges):
                    edges.append((key2, key1))

    return edges

if __name__ == "__main__":

    args = get_args()

    input_file = args.file
    output_file = args.outfile

    fasta_dict = {}

    seq = ""

    with open(input_file) as fh:
        for line in fh:
            if (line[0] == ">") and (seq != ""):
                fasta_dict[header] = seq
                seq = ""
                header = line.strip()[1:]
            elif (line[0] == ">") and (seq == ""):
                header = line.strip()[1:]
            elif line[0] != ">":
                seq += line.strip()
        fasta_dict[header] = seq

    edges = ol_graph(fasta_dict, 3)

    with open(output_file, "w") as fh:
        for pair in edges:
            fh.write(f"{pair[0]} {pair[1]}\n")

