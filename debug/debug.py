#!/usr/bin/env python

# Author: Ian VanGordon
# 12/07/2023

'''
This script tests all of the functions in this repository to ensure they are all working properly
'''

############################################################################################################################################################################################################################################################
# Importing functions
############################################################################################################################################################################################################################################################

import sys
sys.path.append("..")
from bioinfo_toolbox import *


############################################################################################################################################################################################################################################################
# Testing
############################################################################################################################################################################################################################################################

# convert_phred()
assert convert_phred("I") == 40, "wrong phred score for 'I'"
assert convert_phred("C") == 34, "wrong phred score for 'C'"
assert convert_phred("2") == 17, "wrong phred score for '2'"
assert convert_phred("@") == 31, "wrong phred score for '@'"
assert convert_phred("$") == 3, "wrong phred score for '$'"
print("convert_phred() function working properly")

# qual_score()
assert qual_score("IC2@$") == 25
assert qual_score("C2C2") == 25.5
assert qual_score("#I") == 21
assert qual_score("EJ") == 38.5
print("qual_score() function working properly")

#validate_base_seq()
assert validate_base_seq("ATGNac") == True
assert validate_base_seq("AUNag", RNAflag = True) == True
assert validate_base_seq("AGHVBA") == False
assert validate_base_seq("ATVGH", RNAflag = True) == False
print("validate_base_seq() function working properly")

# gc_content()
assert gc_content("CGCGCG") == 1.0
assert gc_content("CGATCGAT") == 0.5
assert gc_content("GATAA") == 0.2
assert gc_content("ATATAT") == 0
print("gc_content() function is working properly")

# calc_median()
assert calc_median([1, 2, 3]) == 2.0
assert calc_median([2, 4, 6, 8]) == 5
assert calc_median([1.5, 2.5, 3.2]) == 2.5
print("calc_median() function is working properly")

# oneline_fasta()
oneline_fasta("bioinfo.test", "bioinfo.test.results")
with open("bioinfo.test.results") as fw:
    numLines = 0
    for lines in fw:
        numLines += 1
assert numLines == 12
print("oneline_fast() function is working properly")

# rev_compliment()
assert rev_compliment("AAA") == "TTT"
assert rev_compliment("AUG") == "CAU"
assert rev_compliment("ACGTACT") == "AGTACGT"
assert rev_compliment("UUACGU") == "ACGUAA"
print("rev_compliment() function is working properly")

# protein_mass_calculator()
assert round(protein_mass_calculator("MAG"), 3) == 259.099
assert round(protein_mass_calculator("DEFNMPQR"), 3) == 1017.434
assert round(protein_mass_calculator("YWVTSSMN"), 3) == 968.406
assert round(protein_mass_calculator("HIKLGHIS"), 3) == 885.518
print("protein_mass_calculator() function is working properly")

# DNA_to_aa
assert dna_to_aa("ATG") == "M"
assert dna_to_aa("ATTCGAGAAATG") == "IREM"
assert dna_to_aa("TGAATTTCG") == "*IS"
assert dna_to_aa("GATCTA") == "DL"
print("dna_to_aa() function is working properly")

# permutation_calc()
assert permutation_calc(3, 3, False) == 6
assert permutation_calc(6, 6, False) == 720
assert permutation_calc(5, 3, False) == 60
assert permutation_calc(10, 5, False) == 30240
print("permutation_calc() function is working properly")