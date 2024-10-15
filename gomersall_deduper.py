#!/usr/bin/env/ python
import  argparse

def get_args(): 
    parser = argparse.ArgumentParser(description="Python program for reference based PCR duplicate removal.")
    parser.add_argument("-f", "--input", help="Absolute file path to sorted SAM file.", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="Absolute file path to deduplicated SAM file.", type=str, required=True)
    parser.add_argument("-u", "--umi", help="File containing list of UMI sequences.", type=str, required=True)
    parser.add_argument("-h", "--help", help="Useful information about the program.", type=str, required=False)
    return parser.parse_args()

INSAM = get_args().input 
OUTSAM = get_args().outfile
UMIS = get_args().umi

# code must be able to run (in a single step) if given a command in the format: ./gomersall_deduper.py -u STL96.txt -f <in.sam> -o <out.sam>

# Pseudocode 

read in list of UMIs 

Open files Input and Output 

Read input file line by line. 

    First read will be included in the output. 

    The fields which are useful for deduplication are: 
        QNAME (col 1) the barcode sequence is at the end of the very end of the string. 
        FLAG (col 2) bitwise flag will show the following info: 
            I am pretty sure we just want the 16th bit to determine if there is coding or non-coding strand.  
        RNAME (col 3) has the chromosome  
        CIGAR string (col 6) got some good info
        POS (col 4) mapping position *adjust based on CIGAR string soft-clipping

Check these in order: 
    chromosome
    Position
        softclipped? before checking position you should look at cigar string. 
        if yes, adjust the position value
    UMI Barcodes match
    strandedness


Before I give this script the input sam file, I will use unix commands to sort the file by Chromosome and then by the position.
    Here is my thought: Since I am only ever writing the FIRST read encountered of all the duplicates. 
    How do I know if I am looking at a new "region" of DNA to compare the read bs while I am reading line-by-line? 

    I know that PCR duplicated will be nearby to each other

    They may be separated by a couple lines which may not match UMI.. 

    For each chromosome and position I will make a dictionary? : Key is (adjusted Position), value is a set of tuples ( UMI barcode, strandedness) seen for that location a.k.a. biological replicate. 
    As soon as I have found a calculated position for a read that differs from the last CALCULATED POSITION then you can throw away that list and write that read to the output file (it is the first of this PCR duplicate you have encountered) 




