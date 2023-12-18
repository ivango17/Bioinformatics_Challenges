#!/usr/bin/env python

# Author: Ian VanGordon
# 12/17/2023

'''
This program takes a string of letters and returns every combination of length r of those letters.
This script solves the Rosalind problem Enumerating k-mers Lexicographically.
'''

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="This program takes a string of space separated characters and returns every combination of r length of those characters in lexicographic order.")
    parser.add_argument("-s", "--string", help="What is the space separated string of letters?", type=str)
    parser.add_argument("-r", "--length", help="What how many characters belong in each combination?", type=str)
    return parser.parse_args()

def combinations(letters, r, new_list = []):
    
    old_len = len(new_list)
    result = []
    for letter in letters:
        if len(new_list) >= len(letters):
            
            for i in range(old_len): 
                result.append(new_list[i] + letter)
        
        else:
            
            result.append(letter)
    if len(result[0]) == r:
        return sorted(result)
    else:
        return combinations(letters, r, result)

if __name__ == "__main__":

    args = get_args()

    letters = args.string.split()
    r = int(args.length)

    sorted_list = combinations(letters, r)

    for comb in sorted_list:
        print(comb)



