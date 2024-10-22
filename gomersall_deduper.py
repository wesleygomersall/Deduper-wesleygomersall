#!/usr/bin/env python
import  argparse
import bioinfo 

# code must be able to run (in a single step) if given a command in the format: ./gomersall_deduper.py -u STL96.txt -f <in.sam> -o <out.sam>

def get_args(): 
    parser = argparse.ArgumentParser(description="Python program for reference based PCR duplicate removal.")
    parser.add_argument("-f", "--input", help="Absolute file path to sorted SAM file.", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="Absolute file path to deduplicated SAM file.", type=str, required=True)
    parser.add_argument("-u", "--umi", help="File containing list of UMI sequences.", type=str, required=True)
    parser.add_argument("-h", "--help", help="Useful information about the program.", type=str, required=False)
    return parser.parse_args()

# Before I give this script the input sam file, I will use unix commands to sort the file by Chromosome and then by the position.

INSAM = get_args().input 
OUTSAM = get_args().outfile
UMIS = get_args().umi

def DNAseqfile_to_set(umilistfile: str, revcomp: bool = False) -> set: 
    """Give this function a <file> of valid DNA sequnces, one sequence per line.
        Returns: a set of those sequences found in <file>. 
        If second parameter is True (default is False), 
        Returns: a set of the reverse complements of sequences in <file>."""
    setofumis = set() 
    with open(UMIS, 'r') as umifile: 
        while True: 
            barcode = umifile.readline() 
            bases = set(barcode) 
            if bases.issubset("AaTtCcGg"): # valid sequence
                # this is where I need to remove whitepace from the line
                barcode.strip()
                if revcomp: # reverse complement the read if revcomp option was set to True (default is False) 
                    barcode = reverse.complement(barcode) 
                setofumis.add(barcode) # add the read to setofumis 
            else: # if this barcode is not a legit sequence then do not add it to the set of barcodes. 
                print(f"The sequence {barcode} is not a valid DNA sequence, it was not added to the set.")
        return(setofumis) 

def line_info(line: str):
    """Docstring: THIS FUNCTION IS NOT COMPLETE"""
    # from line get the following variables: 
        # umi from the very end of the first column entry
        # bflag from column 2 # revcomp = (bflag & 16 == 16) 
        # chrom from column 3
        # pos from column 4
        # cigar from column 6 
    splitupline = line.split()
    chrom = splitupline[2]
    umi = splitupline[0].split(':')[-1] # barcode is the last section of the first column entry, separated by ':'
    
    pos = splitupline[3]

    rev = splitupline[1] & 16 == 16 
    
    cigar = splitupline[5]

    print(f"chromosome: {chrom}, umi: {umi}, rev: {rev}") 
    adjpos = pos

    if rev:
        pass
        # split cigar strand by some letters
        # Add all M, and D lengths (for matches and deletions)
        # Add only the S lengths if they are at the very end of the CIGAR string.
        # Add this total length to the pos to get adjpos.

    else:    # the read is not reverse complemented 
        pass
        # split cigar strand by S
            # if the first element of this split has an M in it then do nothing 
            # else read that value as an int and subtract that int from the pos to get adjpos

    return chrom, adjpos, umi, rev

FWDUMI = DNAseqfile_to_set(UMIS)
REVUMI = DNAseqfile_to_set(UMIS, True)
 
with open(INSAM, 'r') as fin, open(OUTSAM, 'w') as fout: 
    # last_chrom = "Inital Value, Not a chromosome" 
    seenreads = dict()  # Keys: adjusted positions; Value: set of tuples containing the UMI and strands for reads at that position 
    # While Loop for writing the beginning of the file. Look for lines that start with @ and write them. 
    while True: 

        linecontents = fin.readline() # read line
        linesep = linecontents.split() # split line, store the first value of the split
        fout.write(f"{linecontents}\n") 
        if linesep[0] in ['@HD', '@SQ', '@RG']: # if that first value is '@HD', '@SQ', or '@RG'
            continue
        else: # this is the first actual read. Save the position and everything and add them to the dict
            chrom, adjpos, barcode, revstranded = line_info(linecontents) # call the function line_info here
            last_chrom = chrom
            seenreads.setdefault(adjpos, set((barcode, revstranded))) # create a set with the first tuple in it
            break
    
    while True: # now are looping through remaining lines in the file. 
        written = False
        linecontents = fin.readline()

        if linecontents == "": 
            break # this is the end of the file 

        chrom, adjpos, barcode, revstranded = line_info(linecontents) # call the function line_info here

        if revstranded: # check for a valid barcode here 
            if barcode not in REVUMI:
                continue
        else: 
            if barcode not in FWDUMI: 
                continue

        if chrom != last_chrom: # first read on a chromosome must be a new read
            written = True 
            seenreads = dict() # erase dictionary

        else: # on the same chromosome, must check position
            if adjpos not in seenreads.keys(): # new position on chromosome, must be a new read
                written = True 
            else: # may be PCR or biological duplicate
                if (barcode, revstranded) not in seenreads[adjpos]: # this is a biological duplicate
                    written = True
                else: # this read has been seen already and is therefore a PCR duplicate
                    pass 
        if written: 
            fout.write(f"{linecontents}")
            if adjpos in seenreads.keys():
                seenreads.values(adjpos).add((barcode, revstranded)) # add the tuple to the value (set) of the lookup table
            else: 
                seenreads.setdefault(adjpos, (barcode, revstranded)) # add the key and the new set to the lookup table
        last_chrom = chrom



