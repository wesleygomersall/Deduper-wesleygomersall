#!/usr/bin/env python
from collections import defaultdict
import re
import argparse
import bioinfo 

# code must be able to run (in a single step) if given a command in the format: ./gomersall_deduper.py -u STL96.txt -f <in.sam> -o <out.sam>

def get_args(): 
    parser = argparse.ArgumentParser(description="Python program for reference based PCR duplicate removal. Input SAM must be sorted using samtools prior to running this program.")
    parser.add_argument("-f", "--input", help="Absolute file path to sorted SAM file.", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="Absolute file path to deduplicated SAM file.", type=str, required=True)
    parser.add_argument("-u", "--umi", help="File containing list of UMI sequences.", type=str, required=True)
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
            
            if barcode == "": 
                break 

            barcode = barcode.strip().replace('\n', '')
            bases = set(barcode)
            if bases.issubset( ['A', 'a', 'T', 't', 'C', 'c', 'G', 'g']): # valid sequence
                if revcomp: # reverse complement the read if revcomp option was set to True (default is False) 
                    barcode = bioinfo.reverse_complement(barcode) 
                setofumis.add(barcode) # add the read to setofumis 
            else: # if this barcode is not a legit sequence then do not add it to the set of barcodes. 
                print(f"The sequence {barcode} is not a valid DNA sequence, it was not added to the set.")
        return(setofumis) 

def line_info(line: str):
    """This function reads a (read) line from a sam file. 
        It outputs information from that line in this order: 
            Chromosome
            Adjusted position (accounting for CIGAR string)
            Barcode (from the end of the QNAME)
            Bool telling whether the read is aligned to the reverse strand"""
    splitupline = line.split()
    chrom = splitupline[2]
    umi = splitupline[0].split(':')[-1] # barcode is the last section of the first column entry, separated by ':'
    pos = splitupline[3]
    rev = int(splitupline[1]) & 16 == 16 
    cigar = splitupline[5]
    splitcigar = re.findall(r'\d*\D', cigar) 

    if rev:
        sum = 0 
        for i in splitcigar: 
            if 'M' in i: 
                sum += int(i.strip('M'))
            if 'D' in i: # consumes reference (see CIGAR string documentation)
                sum += int(i.strip('D'))
            if 'N' in i: # consumes reference (see CIGAR string documentation) 
                sum += int(i.strip('N'))
        if 'S' in splitcigar[-1]:
            sum += int(splitcigar[-1].strip('S'))
        adjpos = int(pos) + sum - 1

    elif 'S' in splitcigar[0]:    # the read is not reverse complemented 
        adjpos = int(pos) - int(splitcigar[0].strip('S')) 
    else: 
        adjpos = int(pos)

    # print(f"chromosome: {chrom}, pos: {pos}, adjpos: {adjpos}, umi: {umi}, rev: {rev}, cigar: {cigar}, {splitcigar}")
    return chrom, adjpos, umi, rev

UMI = DNAseqfile_to_set(UMIS)
 
with open(INSAM, 'r') as fin, open(OUTSAM, 'w') as fout: 
    seenreads = defaultdict(set)  # Keys: adjusted positions; Value: set of tuples containing the UMI and strands for reads at that position 
    while True: # While Loop for writing the beginning of the file. Look for lines that start with @ and write them. 

        linecontents = fin.readline() # read line
        linesep = linecontents.split() # split line, store the first value of the split
        fout.write(f"{linecontents}") 
        if linesep[0] in ['@HD', '@SQ', '@RG', '@PG', '@CO']: # see part 1.3 SAMv1.pdf sam documentation
            continue
        else: # this is the first actual read. Save the position and everything and add them to the dict
            chrom, adjpos, barcode, revstranded = line_info(linecontents) # call the function line_info here
            last_chrom = chrom
            seenreads.setdefault(adjpos, set()).add((barcode, revstranded)) # create a set with the first tuple in it
            break
   
    while True: # now are looping through remaining lines in the file. 
        written = False
        linecontents = fin.readline()

        if linecontents == "": 
            break # this is the end of the file 

        chrom, adjpos, barcode, revstranded = line_info(linecontents) # call the function line_info here
        val: tuple = (barcode, revstranded)
        
        if barcode not in UMI: 
            last_chrom = chrom 
            continue

        if chrom != last_chrom: # first read on a chromosome must be a new read
            written = True 
            seenreads.clear() # erase dictionary
        
        if chrom == last_chrom: # on the same chromosome, must check position
            if adjpos not in seenreads.keys(): # new position on chromosome, must be a new read
                written = True 
            else: # may be PCR or biological duplicate
                if val in seenreads[adjpos]: # this is a biological duplicate
                    pass
                else: # this read has been seen already and is therefore a PCR duplicate
                    written = True

        if written: 
            fout.write(f"{linecontents}")
            if adjpos in seenreads.keys(): # this is a new biological read for this position
                seenreads[adjpos].add(val)
            else: 
                seenreads.setdefault(adjpos, set()).add(val) # add the key and the new set to the lookup table
        last_chrom = chrom
