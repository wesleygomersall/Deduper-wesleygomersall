#!/usr/bin/env python
import  argparse
import bioinfo 

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

def get_umi(umilistfile: str): 
    # Read through the file with this address
    # File should be text of each UMI, one per line
    # check that each umi is a valid sequence
    # Remove whitespace at the end of each line
    # Add an option of whether you want the reverse complement of each read instead. 
    # return a set of the UMIs
    """ Docstring """
    setofumis = set() 
    with open(UMIS, 'r') as umifile: 
        barcode = umifile.readline() 
        bases = set(barcode) 
        if bases.issubset("AaTtCcGg"): # valid sequence
            print("valid seq") 
            # this is where I need to remove whitepace from the line
            # reverse complement the read if revcomp option was set to True (default is False) 
            # add the read to setofumis 
            return(setofumis) 
        else: 
            return(f"the sequence {barcode} is not a valid DNA sequence")

def line_info(line):

    from line get the following variables: 
        umi from the very end of the first column entry
        bflag from column 2 # revcomp = (bflag & 16 == 16) 
        chrom from column 3
        pos from column 4
        cigar from column 6 

    First determine which strand the read aligned to with rev = bflag & 16 == 16

    if rev:
        reverse_complement(umi) 

        split cigar strand by some letters
        Add all M, and D lengths (for matches and deletions)
        ?Subtract all I lengths? 
        Add only the S lengths if they are at the very end of the CIGAR string.
        Add this total length to the pos to get adjpos.

    else    # the read is not reverse complemented 
        split cigar strand by S
            if the first element of this split has an M in it then do nothing 
            else read that value as an int and subtract that int from the pos to get adjpos

    return chrom, adjpos, umi, rev


# code must be able to run (in a single step) if given a command in the format: ./gomersall_deduper.py -u STL96.txt -f <in.sam> -o <out.sam>

# Pseudocode 

Before I give this script the input sam file, I will use unix commands to sort the file by Chromosome and then by the position.

Open files  INSAM and OUTSAM

    set initial variables: 
        last_chrom 
        dictionary (=dict()) # Keys will be adjusted position (pos - soft clipping). Value will be a set of tuples containing the UMI and strands of the reads at that position. 

    Write the first few lines to the output file automatically untill you start hitting reads..

    begin while loop 

        read line from input

        if the line just read is empty then exit loop
        
        written = False

        Call function: line_info(line)
            returns:umi in correct orientation
                    revstrand is a bool which tells if the read was aligned to the non coding strand.
                    chrom is the chromosome number 
                    adjpos which is the adjusted starting position of the read.

        if chrom != last_chrom # the read must be the first for this chromosome and should be written immediately
            write line to output and go to the next line in the input file.
            written = True 
            erase the dictionary!
        else
            
            # at this point we are still on the same chromosome but now I need to check position along the chromosome 

            if adjpos is not in the keys of dictionary 
                write the line to the output file. 
                written = True 

            else
                
                Get strandedness from bflag:
                    revcomp = bflag & 16 == 16 or something like that.

                if using UMI list then check whether or not umi is an element of that list.
                    if not then there was an error in sequencing the umi and the read should be thrown out
                    probably will record this with a counter too and report it at the end of the program. 

                if the tuple (umi, rev) is an element of the set which is the value assoc. with the real_pos key in dictionary: 
                    This is PCR duplicate. 

                else # this must be a biological duplicate. 
                    write the line to output
                    written = True 

        if written True 
            add the key (real_pos) to dictionary if it is not already there, the value for that key is a set of tuples: (UMI sequence, rev)
                If that key was already there, just add (umi, rev) to the set which is the value for that key. 
            last_chrom = chrom
        

    end loop

