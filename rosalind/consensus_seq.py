#!/usr/bin/env python

# Author: Ian VanGordon
# 12/28/2023

'''
This program takes a FASTA file of sequences of all the same length and returns the consensus sequence by choosing the base per position with the largest count.
This script solves the Rosalind challenge Consensus and Profile.
'''

import argparse
import sys
sys.path.append("..")
from bioinfo_toolbox import oneline_fasta

def get_args():
    parser = argparse.ArgumentParser(description="This program takes a FASTA file of sequences of all the same length and returns the consensus sequence by choosing the base per position with the largest count.")
    parser.add_argument("-f", "--file", help="Where is the file containing the FASTA sequence?", type=str)
    parser.add_argument("-o", "--outputfile", help="What is the destination for the output?", type=str)
    parser.add_argument("-t", "--temp", help="What is the destination for the oneline FASTA file?", type=str)
    return parser.parse_args()

def base_comp(seq_string, base_dict_list):
    '''This funcation takes a DNA sequence as well as a running list of base counts and returns the list with updated base counts per position.'''
    if base_dict_list == []:
        for i in range(len(seq_string)):
            base_dictionary = {"A" : 0, "C" : 0, "G" : 0, "T" : 0}

            match seq_string[i]:
                case "A":
                    base_dictionary["A"] += 1
                case "C":
                    base_dictionary["C"] += 1
                case "G":
                    base_dictionary["G"] += 1
                case "T":
                    base_dictionary["T"] += 1
            base_dict_list.append(base_dictionary)

    else:
        for i in range(len(base_dict_list)):

            match seq_string[i]:
                case "A":
                    base_dict_list[i]["A"] += 1
                case "C":
                    base_dict_list[i]["C"] += 1
                case "G":
                    base_dict_list[i]["G"] += 1
                case "T":
                    base_dict_list[i]["T"] += 1

    return base_dict_list

def consensus_seq(base_dict_list):
    '''This function takes a list of dictionaries respective to a sequence and returns the consensus sequence by choosing bases with the highest count.'''
    con_seq = ""
    for i in range(len(base_dict_list)):
        con_base_max = 0
        for base in base_dict_list[i]:
            if base_dict_list[i][base] > con_base_max:
                con_base = base
                con_base_max = base_dict_list[i][base]
        con_seq += con_base

    return con_seq

if __name__ == "__main__":
    args = get_args()
    
    file = args.file
    temp = args.temp
    output = args.outputfile

    oneline_fasta(file, temp)

    base_dict_list = []
    nt = "ACGT"

    with open(temp) as file:
        for line in file:
            if line[0] != ">":
                seq = line.strip() 
                base_dict_list = base_comp(seq, base_dict_list)

    con_seq = consensus_seq(base_dict_list)

    with open(output, "w") as output:
        output.write(f"{con_seq}")
        for base in nt:
            output.write(f"\n{base}: ")
            for i in range(len(base_dict_list)):
                output.write(f"{base_dict_list[i][base]} ")
