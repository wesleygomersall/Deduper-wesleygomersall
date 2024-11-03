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

Next to complete is implementing the choice between first read or highest quality (or longest for single-end data) 
