# Reference based PCR duplicate removal 

The goal of this assignment is to create a python script which, when given a SAM file of uniquely mapped reads, and a text file containing the known UMIs, removes all PCR duplicates (retain only a single copy of each read). Remember:
- Your Python code can assume a sorted sam file (you *might* need to use `samtools sort` outside of your Python script). Also may assume that all reads are mapped. 
- Account for: 
    - all possible CIGAR strings (including adjusting for soft clipping, etc.)
    - Strand
    - Single-end reads
    - Known UMIs
- Considerations:
    - Millions of reads – avoid loading everything into memory!
    - Be sure to utilize functions appropriately
    - Appropriately comment code and include doc strings
- **CHALLENGE**: In a **separate branch**, implement options for
    - Single-end vs paired-end
    - Known UMIs vs randomers
    - Error correction of known UMIs
    - Choice of duplicate written to file
    
You MUST:
- Write Python 3.12 compatible code
- Include the following argparse options
    - ```-f```, ```--file```: designates absolute file path to sorted sam file
    - ```-o```, ```--outfile```: designates absolute file path to deduplicated sam file
    - ```-u```, ```--umi```: designates file containing the list of UMIs
    - ```-h```, ```--help```: prints a USEFUL help message (see argparse docs)
        - That is, your code must be able to run (in a single step) if given a command in the format:
          ```
          ./<your_last_name>_deduper.py -u STL96.txt -f <in.sam> -o <out.sam>
          ```
- Output the first read encountered if duplicates are found
- Output a properly formatted SAM file
- Name your python script ```<your_last_name>_deduper.py``` and place it in the top level of your repo (that is, not inside a folder)

## Tips from Leslie's introduction: 

After looking at the positions in the alignment output, check UMI barcodes. PCR duplicates should contain the same UMI sequence. 

If two reads are PCR duplicates: 
they will share : 
       Chromosome
       5' start position of read (adjust for soft-clipping) 
       strand (+ or -) 
       UMI 

Be careful with CIGAR String from SAM file. Due to soft clipping the reads may have a different position in the alignment however they are in fact PCR duplicates.

## 2024-10-09

Begin creating algorithm pseudocode for deduplication in main python program file: [gomersall-deduper.py](./gomersall-deduper.py)

from SAMv1 doc: 

Fields of SAM file: 

CIGAR string (6th field) 

| Op | BAM | Description | Consumes query | Consumes reference | 
| ---| ---| --- | --- | --- | 
| M | 0 | alignment match (can be a sequence match or mismatch) | yes | yes | 
| I | 1 | insertion to the reference | yes | no | 
| D | 2 | deletion from the reference | no | yes |  
| N | 3 | skipped region from the reference | no | yes | 
| S | 4 | soft clipping (clipped sequences present in SEQ) | yes | no |
| H | 5 | hard clipping (clipped sequences NOT present in SEQ) | no | no |
| P | 6 | padding (silent deletion from padded reference) | no | no |
| = | 7 | sequence match | yes | yes | 
| X | 8 | sequence mismatch | yes | yes | 

As for the test file, I will edit the provided test.sam from Leslie. 

The required cases for testing: 

1. 3 PCR duplicates which DO NOT need position adjustment from the CIGAR string
	- Chromosome: The same for all three reads 
	- Position: The same for all three reads.  
	- UMI: The same for all three reads.  
	- Strandedness: The same for all three reads. 
	- Written: Only the first read should be written. 

2. 3 PCR duplicates which need position adjustment from the CIGAR string
	- Chromosome: The same for all three reads. 
	- Position: The same ADJUSTED pos for all three reads. The pos should be different, but corrections using CIGAR string should be the same.  
	- UMI: The same UMI from STL96.txt for all three reads. 
	- Strandedness: The same for all three reads. 
	- Written: Only the first read should be written. 

3. Unique read
	- Chromosome: Does not matter
	- Position: unique ADJUSTED position from all other reads in file 
	- UMI: any UMI from list STL96.txt
	- Strandedness: Does not matter
	- Written: Yes

4. 2 Biological duplicates
	- Chromosome: Identical between the two 
	- Position: Identical between the two 
	- UMI: Different between the two, with both coming from the list STL96.txt
	- Strandedness: Same strandedness between the two reads.
	- Written: Yes, both are written.  

5. 2 Reads mapped to the same location (chrom, position, PERHAPS UMI? identical) but the strandedness is different. 
	- Chromosome: Identical between the two 
	- Position: Identical between the two   
	- UMI: Identical between the two, chosen from the list STL96.txt
	- Strandedness: Different 
	- Written: Yes, both are written 

6. The same as #5 except with the umi opposite.
	- Chromosome: Identical between the two 
	- Position: Identical between the two 
	- UMI: Different between the two, both from the list STL96.txt
	- Strandedness: Different between the two 
	- Written: Yes, both are written 

7. The same as #5 but with the reads mapping to the same ADJUSTED position. 
	- Chromosome: Identical between the two reads.  
	- Position: Different position but identical ADJUSTED position between the two (use CIGAR string to specify soft-clipping for one read)  
	- UMI: Identical between the two (chosen from list STL96.txt)  
	- Strandedness: Different between the twoi: different Truth of (bflag & 16 ==16) or something. 
	- Written: Yes, both are written.  

8. Read mapped with incorrect UMI from given list: STL96.txt
	- Chromosome: Does not matter
	- Position: Unique from all other reads 
	- UMI: UMI which does not occur in the file STL96.txt 
	- Strandedness: Does not matter 
	- Written: Not if UMI list is specified in options, otherwise yes.  

9. Two biological duplicates which differ from Pos but are corrected with Cigar string
	- Chromosome: The same for two reads.
	- Position: The pos should be different but ADJUSTED pos should be the same between the two. 
	- UMI: Two separate UMIs from STL96.txt
	- Strandedness: The same strand for both reads.  
	- Written: Yes, both are written. 

10. Two reads which map to the exact same position but different strands.
	- Chromosome: The same for two reads.
	- Position: The pos should be different but ADJUSTED pos should be the same between the two. 
	- UMI: Two separate UMIs from STL96.txt
	- Strandedness: Reads come from different strands. 
	- Written: Yes, both are written. 

11. Four reads mapped to the same position, two are biological duplicates with one PCR duplicate per.  
	- Chromosome: The same for all four reads.
	- Position: The pos can be different but ADJUSTED pos should be the same between the two. 
	- UMI: Two separate UMIs from STL96.txt, twice replicated each.
	- Strandedness: Reads come from the same strand. 
	- Written: Yes, both biological replicates are written, however the PCR duplicates are not. 

12. This is a tricky one. It is actualy two separate mappings. 4 reads: I want to split up the lines of a mapping of biological replicate. R1 is kept, R2 is a whole new mapping, R3 is biological duplicate, R4 is a PCR duplicate. 
	- Chromosome: The same for all four reads.
	- Position: The pos MUST be different but ADJUSTED pos should be the same between the first one  and the last two. 
	- UMI: Different UMIs from the list STL96.txt.
	- Strandedness: Reads come from the same strand. 
	- Written: Yes, both biological replicates are written, as well as the novel sequence between them however the PCR duplicate is not. 


## 2024-10-17

Met with Leslie, Confirmed that for reverse strand reads to adjust the position take the 5 based left most mapped pos, add all Ms, all Ds, and the S at the end then subtract 1 from this value before adding it to the pos.
 Ignore all I becuase these do not consume reference. 


## 2024-11-02

Began challenge. First added an error correcting function for UMI mismatches. This function corrects maximum 2 mismatches.  

## 2024-11-04

Reworked the program to run in two parts: first I loop through the file and determine which lines of the file should be written, then I loop through it once more doing the actual writing.
It is not necessary for writing the first of each unique read in the file, however this allows me to easily implement the options for which read should be kept. 
Next to complete is implementing the choice between first read or highest quality (or longest for single-end data). 

### For paired end data: 

See: http://www.htslib.org/algorithms/duplicate.html

So I should look at the following to determine duplicate: 
1. Read reference and position 
2. Mate reference and position
3. Leftmost (TRUE/FALSE: is this read the lowest aligned position between it and its mate?)
4. Pair orientation: FF, FR, RF, RR 
5. UMI sequences

The 7th column (RNEXT) gives the reference contig of the mate. The 8th column (PNEXT) gives the mapping position of the mate.

Assumptions: 

- Data will be sorted with default `samtools sort` command. 

- All reads in the file are uniquely mapped.

- All paired-end reads share the same QNAME

This is from `samtools sort` documentation: 
> Records with the same name will be ordered according to the values of the READ1 and READ2 flags (see samtools flags). When that flag is also equal, ties are resolved with primary alignments first, then SUPPLEMENTARY, SECONDARY, and finally SUPPLEMENTARY plus SECONDARY. Any remaining ties are reported in the same order as the input data. 

## 2024-11-05

I think I am successfully deduplicating paired-end data at this point, keeping the first of each unique read pair.

Make paired-end test files. Verify with `samtools markdup`?

## 2024-11-08 

I am successfully deduplicating my test files but I now realize that I need to account for the fact that not all read pairs will be right next to each other. I will need to go through my test files and create proper reference between paired reads. 

My error is that there does not seem to be a way to use either samtools collate or samtools sort to sort paired end data such that reads are grouped by name within the same reference.

## 2024-11-09

This is going to be easy now that I understand clearly what constitute paired end reads.

I need to do an enumerate loop through the file. Write all the header lines. When I start writing reads I need to as usual first check the barcodes and correct them if user specifies. 
When I have decided to continue from here then the next thing to do is to look at the split up line.

QNAME of reads are unique except in the case of paired reads here. We can say this about them because they are *uniquely mapped*. 

1. Check one or two barcodes in the QNAME depending on if the reads are paired or not. 
2. If the QNAME is not in the keys of dictionary1 then add it to that dict. 
    - The value for this will look like this: "line_number:adjusted_position:reverse(bool):leftmost(bool):length:meanqual" 
    - If not paired then leftmost is always false? It shouldnt matter I think.
3. If the QNAME is in the keys of both dicts here then maybe return error here.
4. If the QNAME is in the keys of dictionary1 then add it to a second dictionary: dictionary2
    - Same value as above step 2
5. When a chromosome has looped through (this code should be in an if before step 1)
    - Depending on if paired data or not create a new dict of either of the following formats:
    - Single: {{"adjusted_position:reverse(bool):umi": (line_number, score) ... } ... }
    - Paired: {{("R1_adjusted_position:R1_reverse(bool):R1_umi", "R2_adjusted_position:R2_reverse(bool):R2_umi"): (R1line_number, R2line_number, score) ... } ... }
    - Score is a new thing: it is the linescore, lengthscore, or qualityscore. These are calculated to determine which read should be written. 
        - linescore = lines in file - line_number (this will be largest for the first read encountered). 
    - For paired data, use leftmost to determine which dictionary is R1 and which is R2. 
    - Now loop through the entries of this dict. For each key we only want the element of the value set/list which has the highest score. Add the line number(s) which are in these to a list.
    - The list appended to in this step is that list of lines to write. 
    - clear dictionary1 and dictionary2. 

After writing this script it looks like I need to debug the script not picking up the very first read, as well as not writing the final chromosome. It is not immediately clear why I am not writing the first read, but it is clear why I am not picking up the last read. I need to repeat the bit where I write to dictionary 3 then add to the list at the end of the loop. Currently it hits the end of the file and then exits before doing this. 

See the diff between single end data without umi correction, keeping the first read: 
```
25c25
< 03:NS500451:154:HWKTMBGXX:1:11101:25533:1187:ACAGGACA	0	1	74	36	71M	*	0	0	CTTGGTAACTTTCAGAGAATTAGTCACAACTTCTGAAGCAACCACAGTCCATGCAAGTCGACTGGTTTCTC	6AEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEAEEEEEEE<EEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
---
> 01:NS500451:154:HWKTMBGXX:1:11101:24260:1121:ACAGGACA	0	1	74	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
37a38,44
> 19:NS500451:154:HWKTMBGXX:1:11101:10266:1114:TTCGTTCG	40	3	768	36	4S67M	*	0	0	GTGAAACTCGGTGGGATGAGGCGCTCTTTTATATTGAGTTGGGCTGTGCAGGAGTCTTTTCCCACTTCATT	6<EAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEAEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
> 20:NS500451:154:HWKTMBGXX:1:11101:4191:1194:TTCGCCTA	18	3	693	36	10M10D40M11S	*	0	0	GTTCCTCTTGCATTAATGGCTGTGACCCGATGTTCATATTCTAAGCCTTCAATAAGGCCGGTGGAGCGGTA	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
> 21:NS500451:154:HWKTMBGXX:1:11101:53424:50433:CAAGTCGT	0	3	7877	36	2S69M	*	0	0	GTGAAACTCGGTGGGATGAGGCGCTCTTTTATATTGAGTTGGGCTGTGCAGGAGTCTTTTCCCACTTCATT	6<EAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEAEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
> 22:NS500451:154:HWKTMBGXX:1:11101:347:19169:TATAGCGC	0	3	7890	36	15S56M	*	0	0	GTGGGATGAGGCGCTCTTTTATATTGAGTTGGGCTGTGCAGGAGTCTTTTCCCACTTCATTGACGGCGTAG	6<EEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
> 25:NS500451:154:HWKTMBGXX:1:11101:23713:1207:AAGGTACG	0	3	76200105	36	5S66M	*	0	0	GGCTTTGCAGGTGTACTCCCCGCAGTCAACGACCTGGGTTCTCAGGATCTCCAGGCTGGAGACATAGTTCG	/AEEEEEEEEEEEE/EEEEEEEEEEEE/EEEAEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEAEEA<EEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
> 26:NS500451:154:HWKTMBGXX:1:11101:6943:1217:AGCTACCA	0	3	76200105	36	6S65M	*	0	0	CACTTTGCTAGGCTTGCCAATGCCCACAATGTTCTCGGCAGAAACCCTGAACTCATATTCAAGCCCTTCAT	/AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
> 27:NS500451:154:HWKTMBGXX:1:11101:5571:1221:AGGACATG	0	3	76200110	36	10S61M	*	0	0	TTCCAGGTACACAAAAGTCTTCTGAGTAAACAACCTGTACTTTTTGCTACTTCGGATCTGCTTCTTGTCTT	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```

As shown above, I should be writing read 1 rather than 3, and none of the reads on chrom 3 are being written.
