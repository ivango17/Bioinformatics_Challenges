#!/usr/bin/env python

# Author: Ian VanGordon
# 12/08/2023

import argparse
import sys
sys.path.append("..")
from bioinfo_toolbox import DNA_compliments

# def get_args():
#     parser = argparse.ArgumentParser(description="This program generates an output file from a sorted sam file without PCR replicates. This script does not take hard clipping into account.")
#     parser.add_argument("-f", "--file", help="What is the filepath for the txt file to be read?", type=str)
#     parser.add_argument("-o", "--outfile", help="What should this file be called?", type=str)
#     return parser.parse_args()

# args = get_args()

# input_file = args.file
# output_file = args.outfile

pos_length = {}
counter = 0

# with open(input_file) as fh:
#     for line in fh:
#         if line[0] != ">":
#             fw_seq = line.strip()
fw_seq = "TCAATGCATGCGGGTCTATATGCAT"

for i in range(len(fw_seq)):
    for j in range(i, (len(fw_seq) - 1)):
        if fw_seq[i] == DNA_compliments[fw_seq[j]]:
            dist = j - i
            if dist % 2 == 0:
                dist /= 2
                dist = int(dist)
                print(dist)
                for k in range(dist):
                    if fw_seq[i + k] == DNA_compliments[fw_seq[j - k]]:
                        counter += 1
                if counter == dist:
                    pos_length[i + 1] = dist * 2
                    counter = 0
                else:
                    counter = 0

for key in pos_length:
    print(f"{key} {pos_length[key]}")


