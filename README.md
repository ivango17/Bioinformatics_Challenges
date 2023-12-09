# Bioinformatics Challenges
This project is to solve common bioinformatic challenges using Python.

Many of the functions and scripts here are inspired by coursework at [University of Oregon](https://internship.uoregon.edu/bioinformatics) and challenges on [Rosalind](https://rosalind.info/problems/locations/).

## bioinfo_toolbox.py
This script contains important functions pertaining to bioinformatics.

The table below shows the available functions with descriptions:
| Function | Description | Arguments |
| -------- | ----------- | --------- |
| convert_phred() | Takes a phred score and returns the respective Qscore | *letter, val=33* |
| qual_score() | Takes a string of quality scores and returns the average quality score | *phred_score* |
| validate_base_seq() | Takes a sequence and confirms that it is DNA or RNA by returning bool | *seq, RNAflag=False* |
| gc_content() | Takes a DNA or RNA sequence and returns proportion of sequence that is 'G' or 'C' | *seq* |
| calc_median() | Takes a sorted numerical list and returns the median value from that list | *sortedlist* |
| oneline_fasta() | Takes a FASTA file and outputs a FASTA file where every sequence is only one line | *filer, filew='oneline.fa'* |
| rev_compliment() | Takes a sequence of DNA or RNA and returns the reverse compliment sequence | *seq* |
| dna_to_aa() | Takes a sequence of DNA and returns a list of peptides that are encoded from that DNA sequence | *seq* |
| permutation_calc() | Takes n and r and returns number of permutations and an optional list of numeric permutations | *n, r, perm_out=True* |


## Rosalind Challenges
Here is code to solve some of Rosalind's challenges found in the [rosalind folder](./rosalind/).
| Script | Description | Problem Title |
| -------- | ----------- | --------- |
| [point_mutations.py](./rosalind/point_mutations.py) | Calculates the hamming distance between two sequences | [Counting Point Mutations](https://rosalind.info/problems/hamm/) |
| [open_reading_f.py](./rosalind/open_reading_f.py) | Finds all possible polypeptides from a DNA sequence | [Open Reading Frames](https://rosalind.info/problems/orf/) |
| [restriction_sites.py](./rosalind/restriction_sites.py) | Locates restriction sites by finding reverse palindromes in DNA | [Locating Restriction Sites](https://rosalind.info/problems/revp/) |
| [mRNA_poss.py](./rosalind/mRNA_poss.py) | Calculates number of mRNA sequences a protein sequence could have been derived from | [Inferring mRNA from Protein](https://rosalind.info/problems/mrna/) |
| [shared_motif.py](./rosalind/shared_motif.py) | Finds the longest shared motif in a DNA FASTA file | [Finding a Shared Motif](https://rosalind.info/problems/lcsm/) |

## To Do
- [x] Finish bioinfo_toolbox.py
- [ ] Add more Rosalind code
    - [x] Counting Point Mutation
    - [x] Open Reading Frames 
    - [x] Locating Restriction Sites
    - [x] Infering mRNA from Protein
    - [x] Finding a Shared Motif
- [ ] Add to README.md