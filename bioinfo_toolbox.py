#!/usr/bin/env python

# Author: Ian VanGordon
# 12/07/2023
__version__ = "0.1"

'''
This is a set of common bioinformatic functions.

Here are the functions in this program:

1) convert_phred()
2) qual_score()
3) validate_base_seq()
4) gc_content()
5) calc_median()
6) oneline_fasta()
'''

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


############################################################################################################################################################################################################################################################
# Functions
############################################################################################################################################################################################################################################################

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score. Assumes scores are phred+33.'''
    return ord(letter) - 33


def qual_score(phred_score: str) -> float:
    '''Takes a string of qualty scores and returns the average quality score for that string.'''
    total = 0
    for i in range(len(phred_score)):
        total += convert_phred(phred_score[i])
    avg_score = total / len(phred_score)
    return avg_score


def validate_base_seq(seq: str, RNAflag=False) -> bool:
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


def calc_median(sortedlist: list) -> float:
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


def rev_compliment(seq):
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


def protein_mass_calculator(seq):
    '''This function takes a polypeptide sequence and calculates the mass of the protein.'''
    seq = seq.upper()
    protein_mass = 0

    for aa in seq:
        protein_mass += aa_mass[aa]

    return protein_mass


def DNA_to_aa(seq):
    '''This function takes an RNA sequence and translates it to a polypeptide sequence.'''
    seq =  [(seq[i:i + 3]) for i in range(0, len(seq), 3)]
    aa_seq = ""

    for codon in seq:
        aa_seq += codon_to_aa[codon]

    return aa_seq
    



