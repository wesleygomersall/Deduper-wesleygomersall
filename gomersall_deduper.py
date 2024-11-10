#!/usr/bin/env python
from collections import defaultdict
import re
import argparse
import bioinfo 

# Code must be able to run (in a single step) if given a command in the format: ./gomersall_deduper.py -u STL96.txt -f <in.sam> -o <out.sam>
# Before I give this script the input sam file, I will use unix commands to sort the file by Chromosome and then by the position.

def get_args(): 
    parser = argparse.ArgumentParser(description="Python program for reference based PCR duplicate removal. Input SAM must be sorted using samtools prior to running this program. The first of each duplicated read will be written to the output.")
    parser.add_argument("-f", "--input", help="Absolute file path to sorted SAM file.", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="Absolute file path to deduplicated SAM file.", type=str, required=True)
    parser.add_argument("-u", "--umi", help="File containing list of UMI sequences. UMIs will be compared to this list for matching. If no list is given then no error checking of UMI sequences will take place.", type=str, default=None)
    parser.add_argument("-c", "--choice", help="Specify choice of which PCR duplicate to keep. Options are 'first', 'quality', 'longest'. Default is to use first read.", type=str, default='first')
    parser.add_argument("-p", "--paired", help="Specify if data is paired end with 'paired'. Defaults to 'single'.", type=str, default='single')
    parser.add_argument("-e", "--editumi", help="Specify if UMI should be corrected by up to 2 mismatches with 'yes' or 'no'. Does nothing if no UMI file is given. Default is 'no'.", type=str, default='no')
    return parser.parse_args()

def DNAseqfile_to_set(umilistfile: str, revcomp: bool = False) -> set: 
    """Give this function a <file> of valid DNA sequnces, one sequence per line.
        Returns: a set of those sequences found in <file>. 
        If second parameter is True (default is False), 
        Returns: a set of the reverse complements of sequences in <file>."""
    setofumis = set() 
    with open(umilistfile, 'r') as umifile: 
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
    for i in range(tolerance + 1): # first try finding with 0 mismatches, then with 1, and so on, up to the tolerance
        for umi_from_list in validumis: 
            pbases = zip(umi, umi_from_list) # index-wise comparison of 2 sequences 
            # count how many mismatches between bases. If equal to or less than tolerance return that umi 
            if i == sum(pairwise[0] != pairwise[1] for pairwise in pbases): 
                return umi_from_list
    return None

def firstduplicate(file_in:str, umis: list, paired_end: bool = False) -> list: 
    """Returns list of the line numbers for the first of each unique read in SAM file of single-end reads"""
    readstowrite = list()
    seenbarcodes = dict()
    linenum = 0 
    last_chrom = 'unlikely chromosome name'
    countbadumi = 0
    countpcrdup = 0
    countwritten = 0
    readsperchrom = dict() 
    readsthischrom = 0

    with open(INSAM, 'r') as fin: 
        while True: # writing headers
            linecontents = fin.readline()
            linesep = linecontents.split() 
            if linesep[0] in ['@HD', '@SQ', '@RG', '@PG', '@CO']: # see part 1.3 SAMv1.pdf sam documentation
                readstowrite.append(linenum)
                linenum += 1
            else: 
                if paired_end:
                    matecontents = fin.readline()
                break
        
        first_iteration = True  
        while True:

            written = False
            if not first_iteration:
                linecontents = fin.readline()
                if paired_end: # THIS IS NOT CORRECT!
                    matecontens = fin.readline()

            if linecontents == "": 
                readsperchrom.setdefault(last_chrom, readsthischrom)
                break
        
            chrom, adjpos, barcode, revstranded = line_info(linecontents) 

            if paired_end:
                mate_chrom, mate_adjpos, mate_barcode, mate_revstranded = line_info(matecontents) 
                corrected_barcode_1 = nearestumi(barcode.split('^')[0], umis)
                corrected_barcode_2 = nearestumi(barcode.split('^')[1], umis)

                if adjpos < mate_adjpos: # see leftmost samtools markdup documentation
                    # leftmost = is mate ahead of current read 
                    orientation = (revstranded, mate_revstranded) 
                    lineidentifier = f"{adjpos}{mate_adjpos}{orientation}"
                if adjpos > mate_adjpos:
                    orientation = (mate_revstranded, revstranded) 
                    lineidentifier = f"{mate_adjpos}{adjpos}{orientation}"
                
                corrected_barcode = corrected_barcode_1 + '^' + corrected_barcode_2 # this is always the orientation due to same QNAME fields
            
            else:
                corrected_barcode_1 = nearestumi(barcode, umis)
                corrected_barcode_2 = 'Placeholder Barcode'
                corrected_barcode = corrected_barcode_1
                lineidentifier = f"{adjpos}:{revstranded}" # key for dict

            if corrected_barcode_1 == None or corrected_barcode_2 == None:
                last_chrom = chrom
                if paired_end:
                    countbadumi += 2
                else:
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
                    else: 
                        written = True # biological duplicate

            if written:
                countwritten += 1 # should these be 2 for paired end data? 
                readsthischrom += 1
                readstowrite.append(linenum) # add line number directly to list of those to write
                if paired_end: # add mate also 
                    readstowrite.append(linenum + 1) 
                # add read to seenbarcodes to be able to detect PCR duplicates
                if lineidentifier in seenbarcodes.keys():
                    seenbarcodes[lineidentifier].add(corrected_barcode)
                else: 
                    seenbarcodes.setdefault(lineidentifier, set()).add(corrected_barcode) # add key and new set to the table of seen barcodes

            linenum += 1
            if paired_end:
                linenum += 1 # skip over mate read
            first_iteration = False
            last_chrom = chrom

    print(f"Reads removed due to bad UMIs: {countbadumi}")
    print(f"PCR duplicates removed: {countpcrdup}")
    print(f"Reads written to output: {countwritten}")

    print("Chromosome\tNumWritten")
    for key, value in readsperchrom.items():
        print(f"{key}\t{value}")

    return readstowrite 

def line_info2(line_num: int, line: str) -> str:
    """This function reads a (read) line from a sam file. 
        It outputs information from that line in this order: 
        1. line number
        2. umi (this will have both barcodes for PE data)
        3. adjusted position + bool for strandedness
        4. bool if this read is the first in the pair
        5. Match/mismatch length 
        6. length of seq
        7. average qual
            """
    splitupline = line.split()
    chrom = splitupline[2]
    umi = splitupline[0].split(':')[-1] # barcode is the last section of the first column entry, separated by ':'
    pos = splitupline[3]
    rev = int(splitupline[1]) & 16 == 16 
    first = int(splitupline[1]) & 64 == 64 
    cigar = splitupline[5]
    length = len(splitupline[9])
    qual = bioinfo.qual_score(splitupline[10])
    mlength = 0
    splitcigar = re.findall(r'\d*\D', cigar) 
    
    if rev:
        sum = 0 
        for i in splitcigar: 
            if 'M' in i: 
                mlength = int(i.strip('M'))
                sum += mlength
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
    
    if not rev: 
        for i in splitcigar:
            if 'M' in i:
                mlength = int(i.strip('M'))
    
    
    outstring = f"{line_num}:{umi}:{adjpos}{rev}:{first}:{mlength}:{length}:{qual}"
    return outstring

def find_dup(filename: str, umiset: set, correction: bool = False, paired: bool = False, choice: str = 'first') -> list:
    """Returns list of the line numbers for the PCR duplicate read specified by `choice` in SAM file of uniquely aligned reads. 
    If umiset is empty then there is no check for error UMIs. Correction specifies if UMIs are to be corrected by up to two mismatches.
    Deduplicates paired end data if specified by the paired boolean. 
    Choice of keeping 'first', longest match/mismatch 'length', or highest 'quality' PCR duplicate. 
    Keeps the first duplicate in the case of identical length or quality."""
    readstowrite = list()
    dictionary1 = dict()
    dictionary2 = dict()
    dictionary3 = dict()
    last_chrom = 'firstreadoffile'

    totallines = 0
        
    with open(INSAM, 'r') as fin: 
        while True:
            if fin.readline() != '':
                totallines += 1
            else:
                break

    countbadumi = 0
    countpcrdup = 0
    countwritten = 0
    readsperchrom = dict() 
    readsthischrom = 0

    beginning = True

    with open(INSAM, 'r') as fin: 
        for linenum, line in enumerate(fin):

            linesep = line.split()
            if linesep[0] in ['@HD', '@SQ', '@RG', '@PG', '@CO']: 
                readstowrite.append(linenum)
                continue

            if linenum == totallines - 1: # final read of file
                if paired: # check UMI
                    if readid.split(':')[3] == 'True':
                        # this is the first read, take the first UMI in "UMI^UMI" 
                        thisumi = linesep.split('^')[0]
                        thatumi = linesep.split('^')[1]
                    if readid.split(':')[3] == 'False':
                        thisumi = linesep.split('^')[1]
                        thatumi = linesep.split('^')[0]
                else:
                    thisumi = linesep[0].split(':')[-1]
                    thatumi = ''
    
                if umiset != {}:
                    if correction:
                        # correct UMI(s). If None is returned then continue, as this is an invalid umi. inc the bad UMI counter
                        if paired:
                            thisumi = nearestumi(thisumi, umiset)
                            thatumi = nearestumi(thatumi, umiset)
                            if thisumi == None or thatumi == None:
                                countbadumi += 1
                                continue
                        else:
                            thisumi = nearestumi(thisumi, umiset)
                            if thisumi == None:
                                countbadumi += 1
                                continue
                                
                    if not correction:
                        # only check that the UMI(s) are in the list. if not yes for single or yes both for paired then continue, inc the bad umi counter
                        if paired:
                            if (thisumi not in umiset) or (thatumi not in umiset):
                                countbadumi += 1 
                                continue
                        else:
                            if thisumi not in umiset:
                                countbadumi += 1 
                                continue
                # if UMI(s) are good then add the read to dictionary(ies)
                if linesep[0] not in dictionary1.keys():
                    # this is a new QNAME. add this QNAME to dictionary1 as a key with the value in readid
                    dictionary1.setdefault(linesep[0], readid)
                else:
                    if linesep[0] in dictionary2.keys():
                        print("error, this QNAME is already seen twice")
                    # add linesep[0] to dictionary2 as a key with value as readid
                    dictionary2.setdefault(linesep[0], readid)















            if last_chrom != 'firstreadoffile' and linesep[2] != last_chrom or linenum == totallines - 1:
                if paired:
                    # compare dictionary1 and dictionary2. create new dictionary
                    for name in dictionary1.keys():
                        # print(name)
                        
                        # use name to get barcodes. if specified to correct them then do so here
                        bothumis = name.split(':')[-1] # barcode is the last section of the first column entry, separated by ':'

                        if correction:
                            umi1 = nearestumi(bothumis.split('^')[0], umiset)
                            umi2 = nearestumi(bothumis.split('^')[1], umiset)
                        if not correction:
                            umi1 = bothumis.split('^')[0]
                            umi2 = bothumis.split('^')[1]

                        read1 = dictionary1[name].split(':')
                        read2 = dictionary2[name].split(':')
                        
                        # score = pairedscores(read1, read2, choice) # WIP function
                        linescore = max(totallines - int(read1[0]), totallines - int(read2[0]))
                        lengscore = sum(int(read1[4]), int(read2[4]))
                        qualscore = (float(read1[6]) * int(read1[5]) + float(read2[6]) * int(read2[5])) / (int(read1[5]) + int(read2[5]))

                        if read1[3] == 'True': # determine which dictionary contains the first read in the pair
                            newkey: tuple = (read1[2]+umi1, read2[2]+umi2)
                            newvalue: tuple = (int(read1[0]), int(read2[0]), linescore, lengscore, qualscore)
                        if read2[3] == 'True': 
                            newkey: tuple = (read2[2]+umi1, read1[2]+umi2) # the order of UMIs in QNAME is presumed to be in this proper order already
                            newvalue: tuple = (int(read2[0]), int(read1[0]), linescore, lengscore, qualscore)
                            
                        if newkey in dictionary3.keys():
                            # need to add the newvalue tuple to the set which is in the values of dic3
                            dictionary3[newkey].add(newvalue)
                        elif newkey not in dictionary3.keys():
                            # set default and do it so that you can add more tuples to the set later
                            dictionary3.setdefault(newkey, set()).add(newvalue) 
                        
                if not paired:
                    # use only dictionary1 to create new dict
                    for name in dictionary1.keys():
                        umi1 = name.split(':')[-1] # barcode is the last section of the first column entry, separated by ':'

                        if correction:
                            umi1 = nearestumi(umi1, umiset)

                        read1 = dictionary1[name].split(':')
                        
                        linescore = totallines - int(read1[0])
                        lengscore = int(read1[4])
                        qualscore = float(read1[6])

                        newkey: str = read1[2] + umi1
                        newvalue: tuple = (int(read1[0]), 0, linescore, lengscore, qualscore) 
                        # zero in this tuple is placeholder for reading same index as PE above

                        if newkey in dictionary3.keys():
                            # need to add the newvalue tuple to the set which is in the values of dic3
                            dictionary3[newkey].add(newvalue)
                        elif newkey not in dictionary3.keys():
                            # set default and do it so that you can add more tuples to the set later
                            dictionary3.setdefault(newkey, set()).add(newvalue) 


                for setofreads in dictionary3.values():
                    maxseen = 0
                    firstseen = 0
                    for read in setofreads:

                        if paired: 
                            countpcrdup += 2
                        if not paired:
                            countpcrdup += 1
                        if choice == 'first':
                            if int(read[2]) > firstseen:
                                write1 = read[0]
                                write2 = read[1]
                                firstseen = int(read[2])
                            
                        if choice == 'length':
                            if int(read[3]) == maxseen:
                                if int(read[2]) > firstseen: # linescore > firstseen
                                    write1 = read[0]
                                    write2 = read[1]
                                    maxseen = int(read[3])
                                    firstseen = int(read[2])
                            if int(read[3]) > maxseen: # lengscore > maxseen
                                write1 = read[0]
                                write2 = read[1]
                                maxseen = int(read[3])
                                firstseen = int(read[2])

                        if choice == 'quality':
                            if float(read[4]) == maxseen:
                                if int(read[2]) > firstseen: # linescore > firstseen
                                    write1 = read[0]
                                    write2 = read[1]
                                    maxseen = float(read[4])
                                    firstseen = int(read[2])
                            if float(read[4]) > maxseen: # qualscore > maxseen
                                write1 = read[0]
                                write2 = read[1]
                                maxseen = float(read[4])
                                firstseen = int(read[2])

                    readstowrite.append(write1) # add line num for read to write
                    if paired:
                        readstowrite.append(write2) # also add the mate
                        countpcrdup -= 2 
                        countwritten += 2
                        readsthischrom += 2
                    if not paired:
                        countpcrdup -= 1 
                        countwritten += 1
                        readsthischrom += 1
                
                dictionary1 = dict()
                dictionary2 = dict() 
                dictionary3 = dict() 
                    
                # store reads this chrom in a dic and reset counter
                readsperchrom.setdefault(last_chrom, readsthischrom)
                readsthischrom = 0 

            last_chrom = linesep[2] # set before checking UMI in case of bad UMI
            readid = line_info2(linenum, line)

            # check umis
            if paired:
                if readid.split(':')[3] == 'True':
                    # this is the first read, take the first UMI in "UMI^UMI" 
                    thisumi = linesep.split('^')[0]
                    thatumi = linesep.split('^')[1]
                if readid.split(':')[3] == 'False':
                    thisumi = linesep.split('^')[1]
                    thatumi = linesep.split('^')[0]
            else:
                thisumi = linesep[0].split(':')[-1]
                thatumi = ''

            if umiset != {}:
                if correction:
                    # correct UMI(s). If None is returned then continue, as this is an invalid umi. inc the bad UMI counter
                    if paired:
                        thisumi = nearestumi(thisumi, umiset)
                        thatumi = nearestumi(thatumi, umiset)
                        if thisumi == None or thatumi == None:
                            countbadumi += 1
                            continue
                    else:
                        thisumi = nearestumi(thisumi, umiset)
                        if thisumi == None:
                            countbadumi += 1
                            continue
                            
                if not correction:
                    # only check that the UMI(s) are in the list. if not yes for single or yes both for paired then continue, inc the bad umi counter
                    if paired:
                        if (thisumi not in umiset) or (thatumi not in umiset):
                            countbadumi += 1 
                            continue
                    else:
                        if thisumi not in umiset:
                            countbadumi += 1 
                            continue
            # if UMI(s) are good then add the read to dictionary(ies)
            if linesep[0] not in dictionary1.keys():
                # this is a new QNAME. add this QNAME to dictionary1 as a key with the value in readid
                dictionary1.setdefault(linesep[0], readid)
            else:
                if linesep[0] in dictionary2.keys():
                    print("error, this QNAME is already seen twice")
                # add linesep[0] to dictionary2 as a key with value as readid
                dictionary2.setdefault(linesep[0], readid)
            
        readsperchrom.setdefault(last_chrom, readsthischrom)
        readsthischrom = 0 
        print(f"Bad umi count: {countbadumi}")
        print(f"PCR duplicates removed: {countpcrdup}")
        print(f"Reads written: {countwritten}")

        print("Chromosome\tNumWritten")
        for key, value in readsperchrom.items():
            print(f"{key}\t{value}")

        return readstowrite

# 1. Check one or two barcodes in the QNAME depending on if the reads are paired or not. 
# 2. If the QNAME is not in the keys of dictionary1 then add it to that dict. 
    # - The value for this will look like this: "line_number:adjusted_position:reverse(bool):leftmost(bool):length:meanqual" 
    # - If not paired then leftmost is always false? It shouldnt matter I think.
# 3. If the QNAME is in the keys of both dicts here then maybe return error here.
# 4. If the QNAME is in the keys of dictionary1 then add it to a second dictionary: dictionary2
    # - Same value as above step 2
# 5. When a chromosome has looped through (this code should be in an if before step 1)
    # - Depending on if paired data or not create a new dict of either of the following formats:
    # - Single: {{"adjusted_position:reverse(bool):umi": (line_number, score) ... } ... }
    # - Paired: {{("R1_adjusted_position:R1_reverse(bool):R1_umi", "R2_adjusted_position:R2_reverse(bool):R2_umi"): (R1line_number, R2line_number, score) ... } ... }
    # - Score is a new thing: it is the linescore, lengthscore, or qualityscore. These are calculated to determine which read should be written. 
        # - linescore = lines in file - line_number (this will be largest for the first read encountered). 
    # - For paired data, use leftmost to determine which dictionary is R1 and which is R2. 
    # - Now loop through the entries of this dict. For each key we only want the element of the value set/list which has the highest score. Add the line number(s) which are in these to a list.
    # - The list appended to in this step is that list of lines to write. 
    # - clear dictionary1 and dictionary2. 

if __name__ == '__main__':

    INSAM = get_args().input 
    OUTSAM = get_args().outfile
    UMI = DNAseqfile_to_set(get_args().umi)
    if UMI == None:
        UMI = set()
    # PAIRED = True

    CHOICE = get_args().choice
    if CHOICE not in ['first', 'length', 'quality']:
        print("invalid --choice. Please use 'first', 'length', or 'quality'.")
        exit()
    if get_args().paired == 'paired': 
        PAIRED = True
    else:
        PAIRED = False

    if get_args().editumi not in ['yes', 'no']:
        print("invalid --editumi. Please use 'yes' or 'no'.")
        exit()
    CORRECTION = get_args().editumi == 'yes'

    writeme = find_dup(INSAM, UMI, CORRECTION, PAIRED, CHOICE)
        # writeme = firstduplicate(INSAM, UMI, PAIRED)

    # if get_args().choice == 'hi-quality':
        # writeme = hi_qual_duplicate(INSAM, UMI, PAIRED)
        # pass

    # if get_args().choice == 'longest':
        # pass

    with open(INSAM, 'r') as fin, open(OUTSAM, 'w') as fout: 
        for linenum, line in enumerate(fin): 
            if linenum in writeme:
                fout.write(f"{line}")
