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

1. Unique read:
	- Chromosome and position are unique 
2. 2 Biological duplicates:
	- Chromosome and position are identical, but there is a differnt UMI sequence on each one (same strandedness) 

3. 2 Reads mapped to the same location (chrom, position, PERHAPS UMI? identical) but the strandedness is different. 
	- chrom, position, (perhaps) the UMI identical.
	- Need to create a bitwise value for this which gives different value for (bflag & 16 == 16).  VERIFY THIS EXPRESSION BEFORE MAKING THESE
4. The same as above again except with the umi opposite.
	- If did identical UMIs do the same thing but with different UMIs. 
	- Change the chrom and or position from the above test cases also. 
5. 2 Reads mapped to the same ADJUSTED location (chrom, position, PERHAPS UMI? identical) but the strandedness is different. 
	- chrom and (perhaps) the UMI identical.
	- Has different POS but CIGAR string will correct the positions to the same. 
	- Need to create a bitwise value for this which gives different values for (bflag & 16 == 16).  

6. Read mapped with incorrect UMI from given list: STL96.txt
	- If using the UMI list option for the py program this read should not be written to the output.

7. 3 PCR duplicates which DO NOT need position adjustment from the CIGAR string
	- Only the first of these lines should be written. The DNA sequences do not have to be the same but the chrom, position, umi, strandedness must be the same. 
8. 3 PCR duplicates which need position adjustment from the CIGAR string
	- Only the first of these lines should be written. The DNA sequences do not have to be the same but the chrom, umi, strandedness must be the same. 
	- The Positions must be different but will correctly adjust to the same adjusted position using the CIGAR string. 

9. Biological duplicates which differ from Pos but are corrected with Cigar string
	- Same chrom, same adjusted pos
	- different pos, different UMI, different strand. 
	
