#!/usr/bin/env python

# Author: Ian VanGordon
# 12/09/2023

'''
This script was used to generate all possible open reading frames from a DNA string. 
The script solves the challenge "Open Reading Frames" on Rosalind.
'''

import argparse
import sys
sys.path.append("..")
from bioinfo_toolbox import dna_to_aa, rev_compliment

def get_args():
    parser = argparse.ArgumentParser(description="This program generates an output file from a sorted sam file without PCR replicates. This script does not take hard clipping into account.")
    parser.add_argument("-f", "--file", help="What is the filepath for the txt file to be read?", type=str)
    parser.add_argument("-o", "--outfile", help="What should this file be called?", type=str)
    return parser.parse_args()


def reading_frames(DNA_seq):
    '''This function takes a DNA sequence and returns six codon lists.'''
    fw_seq = DNA_seq
    rv_seq = rev_compliment(fw_seq)

    fw_aa = {}
    rv_aa = {}

    for i in range(0,3):
        fw_aa[f"orf{i+1}"] = dna_to_aa(fw_seq[i:])
        rv_aa[f"orf{i+4}"] = dna_to_aa(rv_seq[i:])

    return fw_aa, rv_aa

def protein_string(orf_dict):
    '''Takes a dictionary containing 3 open reading frames for either plus or minus strand and return all possible poly peptides (starting with methionine and ending with a stop codon).'''
    polypeptides = []
    for key in orf_dict:
        for i in range(len(orf_dict[key])):
            if orf_dict[key][i] == "M":
                cur_protein = ""
                for j in range(len(orf_dict[key]) - i):
                    if orf_dict[key][i + j] == "*":
                        frame_pep = (key, cur_protein)
                        polypeptides.append(frame_pep)
                        break
                    else:
                        cur_protein += orf_dict[key][i + j]

    return polypeptides

if __name__ == "__main__":

    uniq_list = []

    args = get_args()

    input_file = args.file
    output_file = args.outfile

    seq = ""

    with open(input_file) as fh:
        for line in fh:
            if line[0] != ">":
                seq += line.strip()

    fw_aa, rv_aa = reading_frames(seq)
    fw_pp = protein_string(fw_aa)
    rv_pp = protein_string(rv_aa)

    # Creating a unique list of polypeptides for outputfile
    for i in range(len(fw_pp)):
        if fw_pp[i][1] not in uniq_list:
            uniq_list.append(fw_pp[i][1])
    for i in range(len(rv_pp)):
        if rv_pp[i][1] not in uniq_list:
            uniq_list.append(rv_pp[i][1])

    with open(output_file, "w") as fh:
        for polypep in uniq_list:
            fh.write(f"{polypep}\n")



