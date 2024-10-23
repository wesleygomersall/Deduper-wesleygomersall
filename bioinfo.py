#!/usr/bin/env python

# Author: Wesley Gomersall wesg@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''
This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework June 2024 - July 2024.
'''

__version__ = "0.3"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

DNA_bases = "ATCGNatcgn"
RNA_bases = "AUCGNaucgn"

def convert_phred(letter: str, phred33: bool = True) -> int:
    """Converts a single character into a phred score. 
    Pass a bool to distinguish between Phred+33 and Phred+64"""    
    if phred33: 
        return ord(letter) - 33
    else:
        return ord(letter) - 64

def qual_score(phred_score: str, phred33: bool = True) -> float:
    """Converts a string of Phred scores into integer scores. 
    Returns the average score for the entire string.
    Pass a bool to distinguish between Phred+33 and Phred+64""" 
    total: int = 0
    for basescore in phred_score: 
        total += convert_phred(basescore, phred33)
    return total/len(phred_score)

def validate_base_seq(seq: str, RNAflag: bool = False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''  
    bases = set(seq)
    return bases.issubset(RNA_bases) if RNAflag else bases.issubset(DNA_bases)

def gc_content(sequence: str):
    '''For string consisting of A,a,C,c,G,g,T,t,N and n, returns the GC content as a value between 0 and 1'''
    assert validate_base_seq(sequence)
    length: int = len(sequence)
    gcount: int = sequence.count('G') + sequence.count('g')
    ccount: int = sequence.count('C') + sequence.count('c')
    acount: int = sequence.count('A') + sequence.count('a')
    tcount: int = sequence.count('T') + sequence.count('t')
    ncount: int = sequence.count('N') + sequence.count('n')
    if length == (gcount + ccount + acount + tcount + ncount):
        return (ccount + gcount)/length
    
def calc_median(lst: list[int | float]) -> float:
    """Returns the median of a sorted list of ints or floats."""
    if (len(lst) % 2) == 0: #even cardinality
        median1 = lst[len(lst)//2 - 1]
        median2 = lst[len(lst)//2]
        median = (median1 + median2)/2 
    if (len(lst) % 2) == 1: #odd cardinality
        median = lst[len(lst)//2]
    return median

def oneline_fasta(fileinput: str, fileoutput: str):
    """Takes string pointing to input fasta file and string naming desired output file.
    Removes newline characters from within the fasta sequences themselves. 
    Output file contains single line headers followed by single line of sequence for each fasta entry. 
    """
    with open(fileinput, 'r') as fin, open(fileoutput, 'w') as fout:
        line = 0
        while True: 
            string = fin.readline().strip()
            if string == "": #if reading empty lines 
                fout.write('\n')
                break
            if string.startswith(">"): #header line
                if line != 0:  # after the first header, 
                    fout.write('\n') # add a new line char before each header to separate from previous entry
                fout.write(f"{string}\n")
            else: #sequence line
                fout.write(string)
            line += 1

def reverse_complement(sequence: str, RNAflag: bool = False) -> str:
    '''Returns the reverse complement of a DNA or an RNA string (default DNA).
    Input and output are both 5' -> 3'.
    '''
    assert validate_base_seq(sequence, RNAflag)
    DNA_comp = {'A':'T', 'a':'t',
                'T':'A', 't':'a',
                'C':'G', 'c':'g',
                'G':'C', 'g':'c',
                'N':'N', 'n':'n'}
    RNA_comp = {'A':'U', 'a':'u',
                'U':'A', 'u':'a',
                'C':'G', 'c':'g',
                'G':'C', 'g':'c',
                'N':'N', 'n':'n'}
    reversecomp: str = ""
    if RNAflag == False:
        for index in range(len(sequence)):
            reversecomp = reversecomp + DNA_comp[sequence[-1*(index + 1)]]
    elif RNAflag:
        for index in range(len(sequence)):
            reversecomp = reversecomp + RNA_comp[sequence[-1*(index + 1)]]
    return reversecomp

def reverse(sequence: str) -> str:
    '''Returns the reverse of a string.
    '''
    reverse: str = ""
    for index in range(len(sequence)):
        reverse = reverse + sequence[-1*(index + 1)]
    return reverse

def check_n(seq: str) -> bool: 
    '''Bool check if a string contains 'n' or 'N'.
    '''
    bases = set(seq.upper())
    if {'N'}.issubset(bases): 
        return True
    else:
        return False
    
if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")

    assert validate_base_seq("AUCGGCuagcnN", True) == True, "validate_base_seq fails for sequence 'AUCGGCuagcnN'"
    assert validate_base_seq("ACTGTGCNTGcaGgtntG") == True, "validate_base_seq fails for sequence 'ACTGTGCNTGcaGgtntG'"
    assert validate_base_seq("ACTGTGCNT_GcagtntG") == False, "validate_base_seq fails for sequence 'ACTGTGCNT_GcagtntG'"
    assert validate_base_seq("ACTGTGCNTxGcagtntG") == False, "validate_base_seq fails for sequence 'ACTGTGCNTxGcagtntG'"
    assert validate_base_seq("ACTGTGCNT GcagtntG") == False, "validate_base_seq fails for sequence 'ACTGTGCNT GcagtntG'"
    print("validate_base_seq function is working properly!")

    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATCGAT") == 0.5
    assert gc_content("GCATTACGnN") == 0.4
    print("Correctly calculated GC content.")

    assert calc_median([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]) == 5.5
    assert calc_median([1, 2, 3, 4, 5, 6, 7, 8, 9]) == 5
    assert calc_median([7]) == 7
    assert calc_median([12, 15.5, 48, 55, 100, 111]) == 51.5
    print("Correctly calculated median of list.")

    assert reverse("ATCG") == "GCTA"
    assert reverse("ATCGAAA") == "AAAGCTA"
    assert reverse("ATCGCGTGCTCTCTCTTCTCT") == "TCTCTTCTCTCTCGTGCGCTA"
    assert reverse_complement("ATCG") == "CGAT"
    assert reverse_complement("ATCGAAA") == "TTTCGAT"
    assert reverse_complement("AUCGAAA", True) == "UUUCGAU" #RNA
    assert reverse_complement("ATCGCGTGCTCTCTCTTCTCT") == "AGAGAAGAGAGAGCACGCGAT"
    print(f"5'-{reverse_complement('AAAAGTCGTCGATTTGCTCGGGGGCTTTC')}-3'")
    print(f"3'-{reverse('AAAAGTCGTCGATTTGCTCGGGGGCTTTC')}-5'")
    print("Successfully reversed and reverse complemented nucleic acid sequences.")

    assert check_n("ACTCGCT") == False
    assert check_n("ACTCGNT") == True
    assert check_n("ACNCGCT") == True
    print("Successfully checked for 'N' bases.")