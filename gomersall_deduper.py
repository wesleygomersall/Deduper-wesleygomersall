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
    The fields which are useful for deduplication are: 
        QNAME (col 1) the barcode sequence is at the end of the very end of the string. 
        FLAG (col 2) bitwise flag will show the following info: 
            I am pretty sure we just want the 16th bit to determine if there is coding or non-coding strand.  
        RNAME (col 3) has the chromosome  
        CIGAR string (col 6) got some good info
        POS (col 4) mapping position *adjust based on CIGAR string soft-clipping

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

        from line also get the following variables: 
            umi from the very end of the first column entry
            bflag from column 2 # revcomp = (bflag & 16 == 16) 
            chrom from column 3
            pos from column 4
            cigar from column 6 

        if chrom != last_chrom # the read must be the first for this chromosome and should be written immediately
            write line to output and go to the next line in the input file.
            written = True 
            erase the dictionary!
        else
            
            # at this point we are still on the same chromosome but now I need to check position along the chromosome 
            First determine whether the read I am looking at is soft clipped. use the cigar string.
                Something like split the cigar string by `S` then hit the first element of that split with a substing by `M`. do nothing if true
                if not true then I should read the first element as an integer and use the value to correct the pos (subtract this int from pos)

            Determine what the REAL position of the read is by subtracting the correction from pos.

            if real_pos is not in the keys of dictionary 
                write the line to the output file. 
                written = True 

            else
                
                Get strandedness from bflag:
                    revcomp = bflag & 16 == 16 or something like that.

                if the read is reverse complemented then I need to Rev complement the umi sequence. 

                if using UMI list then check whether or not umi is an element of that list.
                    if not then there was an error in sequencing the umi and the read should be thrown out
                    probably will record this with a counter too and report it at the end of the program. 

                if the tuple (umi, strandedness) is an element of the set which is the value assoc. with the real_pos key in dictionary: 
                    This is PCR duplicate. 

                else # this must be a biological duplicate. 
                    write the line to output
                    written = True 

        if written True 
            add the key (real_pos) to dictionary if it is not already there, the value for that key is a set of tuples: (UMI sequence, strandedness)
                If that key was already there, just add (umi, strandedness) to the set which is the value for that key. 
            last_chrom = chrom
        

    end loop

