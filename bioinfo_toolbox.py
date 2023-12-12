#!/usr/bin/env python

# Author: Ian VanGordon
# 12/07/2023
__version__ = "0.3"

'''
This is a set of common bioinformatic functions.

Here are the functions in this program:

1) convert_phred()
2) qual_score()
3) validate_base_seq()
4) gc_content()
5) calc_median()
6) oneline_fasta()
7) permutation_calc()
8) transition_transversion()
'''

from math import factorial
from random import sample

############################################################################################################################################################################################################################################################
# Constants
############################################################################################################################################################################################################################################################

DNA_bases = set("AGCTNagctn")
RNA_bases = set("AGCUNagcun")

DNA_compliments = {"A":"T", "T": "A", "G":"C", "C":"G"}
RNA_compliments = {"A":"U", "U": "A", "G":"C", "C":"G"}

aa_mass = {
'A': 71.03711, 'G': 57.02146, 'M': 131.04049, 'S': 87.03203,
'C': 103.00919, 'H': 137.05891, 'N': 114.04293, 'T': 101.04768,
'D': 115.02694, 'I': 113.08406, 'P': 97.05276, 'V': 99.06841,
'E': 129.04259, 'K': 128.09496, 'Q': 128.05858, 'W': 186.07931,
'F': 147.06841, 'L': 113.08406, 'R': 156.10111, 'Y': 163.06333}

codon_to_aa = {
    'TCA': 'S','TCC': 'S','TCG': 'S','TCT': 'S', 'TTC': 'F','TTT': 'F',
    'TTA': 'L','TTG': 'L','TAC': 'Y','TAT': 'Y','TAA': '*','TAG': '*',
    'TGC': 'C','TGT': 'C','TGA': '*','TGG': 'W','CTA': 'L','CTC': 'L',
    'CTG': 'L','CTT': 'L','CCA': 'P','CCC': 'P','CCG': 'P','CCT': 'P',
    'CAC': 'H','CAT': 'H','CAA': 'Q','CAG': 'Q','CGA': 'R','CGC': 'R',
    'CGG': 'R','CGT': 'R','ATA': 'I','ATC': 'I','ATT': 'I','ATG': 'M',
    'ACA': 'T','ACC': 'T','ACG': 'T','ACT': 'T','AAC': 'N','AAT': 'N',
    'AAA': 'K','AAG': 'K','AGC': 'S','AGT': 'S','AGA': 'R','AGG': 'R',
    'GTA': 'V','GTC': 'V','GTG': 'V','GTT': 'V','GCA': 'A','GCC': 'A',
    'GCG': 'A','GCT': 'A','GAC': 'D','GAT': 'D','GAA': 'E','GAG': 'E',
    'GGA': 'G','GGC': 'G','GGG': 'G','GGT': 'G'}

aa_numcodon = {
    'A': 4,'C': 2,'D': 2,'E': 2,
    'F': 2,'G': 4,'I': 3,'H': 2,
    'K': 2,'L': 6,'M': 1,'N': 2,
    'P': 4,'Q': 2,'R': 6,'S': 6,
    'T': 4,'V': 4,'W': 1,'Y': 2,
    '*': 3,}

transversions = {
    'A' : set("CT"), 'G' : set("CT"), 'T' : set("AG"), 'C' : set("AG")}


############################################################################################################################################################################################################################################################
# Functions
############################################################################################################################################################################################################################################################

def convert_phred(letter: str, val: int = 33):
    '''Converts a single character into a phred score. Assumes scores are phred+33.'''
    return ord(letter) - val


def qual_score(phred_score: str):
    '''Takes a string of qualty scores and returns the average quality score for that string.'''
    total = 0
    for i in range(len(phred_score)):
        total += convert_phred(phred_score[i])
    avg_score = total / len(phred_score)
    return avg_score


def validate_base_seq(seq: str, RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case
    insensitive.'''
    seqset = set(seq)
    return seqset <= (RNA_bases if RNAflag else DNA_bases)


def gc_content(seq: str):
    '''Validates that seq is DNA/RNA, then returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(seq) == True
    seq = seq.upper()
    numGC = 0
    for i in range(len(seq)):
        if seq[i] == "G" or seq[i] == "C":
            numGC += 1
    GCcontent = numGC / len(seq)
    return GCcontent


def calc_median(sortedlist: list):
    '''This function finds the median index of a sorted list and returns the value at that index.'''
    if len(sortedlist) % 2 != 0:
        medianIndex = len(sortedlist) // 2
        return sortedlist[medianIndex]
    else:
        i1 = len(sortedlist) // 2
        i2 = i1 - 1
        medianEvens = (sortedlist[i1] + sortedlist[i2]) / 2
        return medianEvens


def oneline_fasta(filer: str, filew: str = "oneline.fa"):
    '''Takes a FASTA file with multiple lines of sequence per read and reduces it to one line of sequence per read'''
    with(open(filer) as fr, open(filew, "w") as fw):
        seq = ""
        header = ""
        for line in fr:
            line = line.strip("\n")
            if line[0] == ">" and header == "":
                header = line
            elif line[0] == ">":
                fw.write(f"{header}\n{seq}\n")
                seq = ""
                header = line
            else:
                seq += line
        fw.write(f"{header}\n{seq}\n")


def rev_compliment(seq: str):
    '''This function takes a sequence of DNA or RNA and reverse compliments the sequence'''
    seq = seq.upper()[::-1]
    new_seq = ""
    if "U" in seq:
        comp = RNA_compliments
    else:
        comp = DNA_compliments
    for bp in seq:
        new_seq += comp[bp]
    return new_seq


def protein_mass_calculator(seq: str):
    '''This function takes a polypeptide sequence and calculates the mass of the protein.'''
    seq = seq.upper()
    protein_mass = 0
    for aa in seq:
        protein_mass += aa_mass[aa]
    return protein_mass


def dna_to_aa(seq: str):
    '''This function takes an RNA sequence and translates it to a polypeptide sequence.'''
    seq =  [(seq[i:i + 3]) for i in range(0, len(seq), 3)]
    aa_seq = ""
    for codon in seq:
        if len(codon) % 3 == 0:
            aa_seq += codon_to_aa[codon]
    return aa_seq
    

def permutation_calc(n: int, r: int, perm_out: bool = True):
    '''This function returns the number of permutations given n and r. Outputs and optional list of all possible permutations.'''
    num_permutations = int(factorial(n) / factorial(n-r))
    if perm_out:
        perm_list = []
        cur_list = []
        population = [(i + 1) for i in range(n)]
        while len(perm_list) != num_permutations:
            cur_list = sample(population, k = r)
            if cur_list not in perm_list:
                perm_list.append(cur_list)
        return num_permutations, perm_list

    else:
        return num_permutations


def transition_transverion(seq1, seq2):
    '''This function takes two DNA sequences and returns the transition to transversion ratio.'''
    transitions_count = 0
    transversions_count = 0

    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            if seq2[i] in transversions[seq1[i]]:
                transversions_count += 1
            else:
                transitions_count += 1

    return transitions_count / transversions_count


def fibonacci_seq(n, k, count=1, f1=0, f2=1):
    '''The following function returns a fibonacci number given n and k(mult per instance).'''
    if count == n:
        return f2
    else:
        temp = f2
        f2 += f1 * k
        f1 = temp
        return fibonacci_seq(n, k, count + 1, f1, f2)

print(fibonacci_seq(5, 3))

