# Deduper

## Part 1
Use this repo template to create your own Deduper repo - you should do all your work in your own repository. Please name it `Deduper-<github-user-name>`.

Write up a strategy for writing a Reference Based PCR Duplicate Removal tool. That is, given a sorted sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read). Develop a strategy that avoids loading everything into memory. You should not write any code for this portion of the assignment. Be sure to:
- Define the problem
- Write examples:
    - Include a properly formated sorted input sam file
    - Include a properly formated expected output sam file
- Develop your algorithm using pseudocode
- Determine high level functions
    - Description
    - Function headers
    - Test examples (for individual functions)
    - Return statement
    
For this portion of the assignment, you should design your algorithm for single-end data, with 96 UMIs. UMI information will be in the QNAME, like so: ```NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT```. Discard any UMIs with errors (or think about how you might error correct, if you're feeling ambitious).

## Part 2
An important part of writing code is reviewing code - both your own and other's. In this portion of the assignment, you will be assigned 3 students' pseudocode algorithms to review. Be sure to evaluate the following points:
- Does the proposed algorithm make sense to you? Can you follow the logic?
- Does the algorithm do everything it's supposed to do? (see part 1)
- Are proposed functions reasonable? Are they "standalone" pieces of code?

You can find your assigned reviewees on Canvas. You can find your fellow students' repositories at 
```
github.com/<user>/Deduper-<github-user-name>
```
Be sure to leave comments on their repositories by creating issues or by commenting on the pull request.

## Part 3
Write your deduper function!

Given a SAM file of uniquely mapped reads, and a text file containing the known UMIs, remove all PCR duplicates (retain only a single copy of each read). Remember:
- Your Python code can assume a sorted sam file (you *might* need to use `samtools sort` outside of your Python script)
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


This is the challenge branch for barcode sensitive, reference-based deduplication of uniquely mapped reads in SAM format. 

The script `gomersall_deduper.py` handles both paired-end and single-end data. User has the choice of specifying which PCR duplicate should be written to output sam file. User also has the choice of whether or not to specify a list of valid UMI sequences to write. These sequences will be compared to the barcodes which are found at the end of the QNAME for each read. Paired end barcodes must be in the format `AAAAAA^TTTTTT`. If no list is given then there is no checking of UMI barcodes for validity. User may also specify whether UMI sequences may be corrected by up to 2 mismatches from those supplied by the user. 

Instructions for running this script: 

`./gomersall_deduper.py -f <input.sam> -o <output.sam> -u <UMIlist.txt> -c 'choice' -p 'paired' -e 'no'`

Explanation of options: 
- `-f`: Input SAM file. Must be sorted with `samtools sort` before being passed into this program. 
- `-o`: Output SAM file. Will have PCR duplicates removed.
- `-u`: Text file containing sequences of valid UMIs, one per line. If there is no list passed, then UMIs in input SAM will not be checked for validity. 
- `-c`: Choice of 'first', 'quality', or 'length'. This specifies which PCR duplicate should be written. The 'quality' and 'length' options will specify the highest average quality read or the longest match/mismatch (in CIGAR string) PCR duplicate to be written, respectively. In case of the same quality or length, the first of these will be written.  
- `-p`: Specify 'paired' or 'single' end data in the file. 
- `-e`: Specify with 'yes' or 'no' whether UMIs should be corrected by up to 2 mismatches, to the closest sequence in the list of specified UMIs. 

