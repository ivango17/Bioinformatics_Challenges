#!/usr/bin/env python

# Author: Ian VanGordon
# 12/12/2023

'''
This script takes n being number of iterations and k multiplier for each iteration.
This program solves the Rosalind challenge Rabbits and Reoccurance Relation.
'''

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="This script takes n being number of iterations and k multiplier for each iteration.")
    parser.add_argument("-n", "--generations", help="How many months?", type=str)
    parser.add_argument("-k", "--lit_count", help="How many offspring per generation?", type=str)
    return parser.parse_args()

def fibonacci_seq(n, k, count=1, f1=0, f2=1):
    '''The following function returns a fibonacci number given n and k(mult per instance).'''
    if count == n:
        return f2
    else:
        temp = f2
        f2 += f1 * k
        f1 = temp
        return fibonacci_seq(n, k, count + 1, f1, f2)

if __name__ == "__main__":
    args = get_args()

    n = int(args.generations)
    k = int(args.lit_count)

    print(fibonacci_seq(n, k))
    