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



