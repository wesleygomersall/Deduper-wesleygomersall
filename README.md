# Deduper challenge 

This is the challenge branch for barcode sensitive, reference-based deduplication of uniquely mapped reads in SAM format. 

The script `gomersall_deduper.py` handles paired-end or single-end data. User has the choice of specifying which PCR duplicate should be written to output sam file. User also has the choice of whether or not to specify a list of valid UMI sequences to write. These sequences will be compared to the barcodes which are found at the end of the QNAME for each read. Paired end barcodes must be in the format `AAAAAA:TTTTTT` or `AAAAAA^TTTTTT`.

## Running this script: 

`./gomersall_deduper.py -f <input.sam> -o <output.sam> -u <UMIlist.txt> -c 'quality' -p 'paired' -e 'no'`

Explanation of options: 
- `-f`: Input SAM file. Must be sorted with `samtools sort` before being passed into this program. 
- `-o`: Output SAM file. Will have PCR duplicates removed.
- `-u`: Text file containing sequences of valid UMIs, one per line. If there is no list passed, then UMIs in input SAM will not be checked for validity. 
- `-c`: Choice of 'first', 'quality', or 'length'. This specifies which PCR duplicate should be written. The 'quality' and 'length' options will specify the highest average quality read or the longest match/mismatch (in CIGAR string) PCR duplicate to be written, respectively. In case of the same quality or length, the first of these will be written.  
- `-p`: Specify 'paired' or 'single' end data in the file. 
- `-e`: Specify with 'yes' or 'no' whether UMIs should be corrected by up to 2 mismatches, to the closest sequence in the list of specified UMIs. Should be 'no' or not specified when not using option `-u`.
