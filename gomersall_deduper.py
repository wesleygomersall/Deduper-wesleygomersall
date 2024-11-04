#!/usr/bin/env python
from collections import defaultdict
import re
import argparse
import bioinfo 

# code must be able to run (in a single step) if given a command in the format: ./gomersall_deduper.py -u STL96.txt -f <in.sam> -o <out.sam>

def get_args(): 
    parser = argparse.ArgumentParser(description="Python program for reference based PCR duplicate removal. Input SAM must be sorted using samtools prior to running this program. The first of each duplicated read will be written to the output.")
    parser.add_argument("-f", "--input", help="Absolute file path to sorted SAM file.", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="Absolute file path to deduplicated SAM file.", type=str, required=True)
    parser.add_argument("-u", "--umi", help="File containing list of UMI sequences. UMIs will be compared to this list for matching and/or error correction up to two base mismatches.", type=str, required=True)
    parser.add_argument("-c", "--choice", help="Specify choice of which PCR duplicate to keep. Options are 'first', 'hi-quality', 'longest'. Default is to use first read. Longest is only applicable to single-end data.", type=str, default='first')
    return parser.parse_args()

# Before I give this script the input sam file, I will use unix commands to sort the file by Chromosome and then by the position.

INSAM = get_args().input 
OUTSAM = get_args().outfile
UMIS = get_args().umi
KEEP = get_args().choice

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

def nearestumi(umi: str, validumis: list, tolerance: int = 2):
    """For some UMI, checks through the list of valid sequences and returns the nearest match (default up to 2 mismatches). 
    If there is no sequence within tolerance, returns 'None'.""" 
    for umi_from_list in validumis: 
        pbases = zip(umi, umi_from_list) # index-wise comparison of 2 sequences 
        # count how many mismatches between bases. If equal to or less than tolerance return that umi 
        if tolerance <= sum(seen == listed for seen, listed in pbases): 
            return umi_from_list
    return None

def filter_criteria(line: str, criteria: str) -> int: 
    """ For criteria choice between 'longest', or 'hi-quality', will return either the appropriate information used to compare reads to keep.
    If 'first' returns none, as this information cannot be found given only a single line. """
    if criteria == 'first': 
        return None
    elif criteria == 'hi-quality': 
        qualstring = line.split()[-1]
        result = qual_score(qualstring) # average quality of read
    elif criteria == 'longest':
        # WIP
        # calculate how many bases were mapped for the read (look at CIGAR string for this one) 
        result = "length"
    return result

UMI = DNAseqfile_to_set(UMIS)


def hi_qual_duplicate(file_in: str, umis: list) -> list:
    readcompare = dict() # dict = {{uniqueread info, info about read(depends on option for --choice}, ... }  
    readlinenum = dict() # dict = {{uniqueread info, line number of read to keep}, ... }  
    readstowrite = list() # store values and clear readlinenum 
    linenum: int = 0 
    criteria = 'hi-quality'

    with open(INSAM, 'r') as fin:
    
    while True: # writing headers
        linecontents = fin.readline()
        linesep = linecontents.split() 
        if linesep[0] in ['@HD', '@SQ', '@RG', '@PG', '@CO']: # see part 1.3 SAMv1.pdf sam documentation
            readstowrite.append(linenum)
            linenum += 1
        else: 
            linenum += 1 
            break

    # at this point we are looking at the first read of the file. 
    # look at linesep and store into readcompare

    chrom, adjpos, barcode, revstranded = line_info(linecontents) 
    readidentifier = f"{adjpos}:{barcode}:{revstranded}" 

    readstat = filter_criteria(linecontents, criteria)
    
    readcompare.setdefault(readidentifier, readstat) # create a set with the first tuple in it
    readstowrite.setdefault(readidentifier, linenum) 
    # store into readlinenum


def firstduplicate(file_in:str, umis) -> list: 
    readstowrite = list()
    seenbarcodes = dict()
    linenum = 0 
    last_chrom = 'unlikely chromosome name'
    countbadumi = 0
    countpcrdup = 0
    countwritten = 0
    readsthischrom = 0

    with open(INSAM, 'r') as fin: 
        while True: # writing headers
            linecontents = fin.readline()
            linesep = linecontents.split() 
            if linesep[0] in ['@HD', '@SQ', '@RG', '@PG', '@CO']: # see part 1.3 SAMv1.pdf sam documentation
                readstowrite.append(linenum)
                linenum += 1
            else: 
                linenum += 1 
                break

        first_iteration = True  
        while True: 
            written = False
            linecontents = fin.readline()

            if linecontents == "": 
                readsperchrom.setdefault(last_chrom, readsthischrom)
                break
        
            chrom, adjpos, barcode, revstranded = line_info(linecontents) # call the function line_info here
            lineidentifier = f"{adjpos}:{revstranded}" # key for dict
            
            corrected_barcode = nearestumi(barcode, UMI)
            if corrected_barcode == None:
                last_chrom = chrom 
                countbadumi += 1
                continue
    
            if chrom != last_chrom: 
                written = True # first read on the chromosome must be the first of its kind in the file
                if not first_iteration:
                    readsperchrom.setdefault(last_chrom, readsthischrom)
                    readsthischrom = 0
                seenbarcodes.clear()
    
            if chrom == last_chrom: 
                if lineidentifier not in seenbarcodes.keys():
                    written = True
                else: 
                    if corrected_barcode in seenbarcodes[lineidentifier]:
                        countpcrdup += 1
                        pass # barcode already seen at this position => PCR duplicate
                    else: 
                        written = True # biological duplicate
            if written:
                countwritten += 1
                readsthischrom += 1
                # add line number directly to list of those to write
            linenum += 1
            first_iteration = False
            last_chrom = chrom

#         if written: 
#             countwritten += 1
#             readsthischrom += 1
#             fout.write(f"{linecontents}")
#             if adjpos in seenreads.keys(): # this is a new biological read for this position
#                 seenreads[adjpos].add(val)
#             else: 
#                 seenreads.setdefault(adjpos, set()).add(val) # add the key and the new set to the lookup table
#         last_chrom = chrom
#         first_iteration = False

    print(f"Reads removed due to bad UMIs: {countbadumi}")
    print(f"PCR duplicates removed: {countpcrdup}")
    print(f"Reads written to output: {countwritten}")

    print("Chromosome\tNumWritten")
    for key, value in readsperchrom.items():
        print(f"{key}\t{value}")










with open(INSAM, 'r') as fin, open(OUTSAM, 'w') as fout: 
    for linenum, line in enumerate(fin): 
        if linenum in readstowrite: 
            fout.write(f"{line}")



# with open(INSAM, 'r') as fin, open(OUTSAM, 'w') as fout: 
#     seenreads = defaultdict(set)  # Keys: adjusted positions; Value: set of tuples containing the UMI and strands for reads at that position 
# 
#     while True: # writing headers 
#         linecontents = fin.readline() 
#         linesep = linecontents.split() 
#         fout.write(f"{linecontents}") 
#         if linesep[0] in ['@HD', '@SQ', '@RG', '@PG', '@CO']: # see part 1.3 SAMv1.pdf sam documentation
#             continue
#         else: # the first actual read. Save everything and add to the dict
#             chrom, adjpos, barcode, revstranded = line_info(linecontents) 
#             last_chrom = chrom
#             seenreads.setdefault(adjpos, set()).add((barcode, revstranded)) # create a set with the first tuple in it
#             break
#     
#     first_iteration = True
#     while True: # rest of file. 
#         written = False
#         linecontents = fin.readline()
# 
#         if linecontents == "": 
#             readsperchrom.setdefault(last_chrom, readsthischrom)
#             break 
# 
#         chrom, adjpos, barcode, revstranded = line_info(linecontents) # call the function line_info here
#         val: tuple = (barcode, revstranded)
#         
#         corrected_barcode = nearestumi(barcode, UMI)
#         if corrected_barcode == None:
#             last_chrom = chrom 
#             countbadumi += 1
#             continue
# 
#         if chrom != last_chrom: 
#             written = True # first read on a chromosome must be a new read
#             if not first_iteration: 
#                 readsperchrom.setdefault(last_chrom, readsthischrom)
#                 readsthischrom = 0
#             seenreads.clear() 
#         
#         if chrom == last_chrom: 
#             if adjpos not in seenreads.keys(): 
#                 written = True # new position on chromosome, must be a new read
#             else: 
#                 if val in seenreads[adjpos]: 
#                     countpcrdup += 1
#                     pass # read has been seen already therefore is a PCR duplicate
#                 else: 
#                     written = True # biological duplicate
# 
#         if written: 
#             countwritten += 1
#             readsthischrom += 1
#             fout.write(f"{linecontents}")
#             if adjpos in seenreads.keys(): # this is a new biological read for this position
#                 seenreads[adjpos].add(val)
#             else: 
#                 seenreads.setdefault(adjpos, set()).add(val) # add the key and the new set to the lookup table
#         last_chrom = chrom
#         first_iteration = False
 
# write the output file after figuring out which reads should be kept. 
