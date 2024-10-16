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

    For each chromosome and ADJUSTED position I will make a dictionary? : 
        Key is tuple (UMI barcode, strand from bitwise flag ) (check if in list if using UMI file), value is a count of all the duplicates seen for that location a.k.a. biological replicate. 
    As soon as I have found a calculated position for a read that differs from the last CALCULATED POSITION then you can throw away that list and write that read to the output file (it is the first of this PCR duplicate you have encountered. 




Open files  
    set initial variables: 
        last_chrom 
        last_position,
        dictionary. (=dict())
    begin loop 
        read line from input

        if the line just read is empty then exit loop
        
        written = False

        from line also get the following variables: 
            umi from the very end of the first column entry
            bflag from column 2
            chrom from column 3
            cigar from column 6 
            pos from column 4

        if chrom != last_chrom # the read must be the first for this chromosome and should be written immediately
            write line to output and go to the next line in the input file.
            written = True 
            erase the dictionary!
        else
            
            # at this point we are still on the same chromosome but now I need to check position along the chromosome 
            First determine whether the read I am looking at is soft clipped. This uses the cigar string.
                Something like split the cigar string by S then hit the first element of that split with a substing M. Pass if true
                    if not true then I should read the first element as an integer and use the value to correct the pos (subtract this int from pos)
            Determine what the REAL position of the read is by subtracting the correction from pos.

            if real_pos != last_position # this is a new location on the chromosome, and the first occurance should be kept.
                write the line to the output file. 
                written = True 
                erase the dictionary!
            else
                
                Get strandedness from bflag:
                    revcomp = bflag & 16 == 16 or something like that.

                if the read is reverse complemented then I need to Rev complement the umi sequence. 

                If (umi, revcomp) is in keys of dictionary 
                    this is a pcr duplicate
                    inc the value by one for this key
                else
                    write the line to output
                    written = True 

        if written True 
            add the key to dictionary

            last_chrom = chrom
            last_position = real_pos
        

                



    end loop

