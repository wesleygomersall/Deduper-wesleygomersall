
20 or 21 reads in output file depending on read #16

Sam header:

```
@HD	VN:1.0	SO:unsorted
@PG	ID:GSNAP	PN:gsnap	VN:2017-10-12	CL:gsnap.avx2 --gunzip -t 26 -A sam -m 5 -d mm10_chroms -D /projects/bgmp/coonrod/mmu/INTEL -s /projects/bgmp/coonrod/mmu/INTEL/mm10_chroms/mm10_chroms.maps/Mus_musculus.GRCm38.89.splicesites.iit --split-output=/projects/bgmp/coonrod/deduper/gsnap//Datset1 /projects/bgmp/coonrod/deduper//Dataset1.fastq_dups.gz
@SQ	SN:1	LN:195471971
@SQ	SN:2	LN:182113224
@SQ	SN:3	LN:160039680
@SQ	SN:4	LN:156508116
@SQ	SN:5	LN:151834684
@SQ	SN:6	LN:149736546
@SQ	SN:7	LN:145441459
@SQ	SN:8	LN:129401213
@SQ	SN:9	LN:124595110
@SQ	SN:10	LN:130694993
@SQ	SN:11	LN:122082543
@SQ	SN:12	LN:120129022
@SQ	SN:13	LN:120421639
@SQ	SN:14	LN:124902244
@SQ	SN:15	LN:104043685
@SQ	SN:16	LN:98207768
@SQ	SN:17	LN:94987271
@SQ	SN:18	LN:90702639
@SQ	SN:19	LN:61431566
@SQ	SN:X	LN:171031299
@SQ	SN:Y	LN:91744698
@SQ	SN:MT	LN:16299
```

# Reads

1-3:
3 PCR duplicates which DO NOT need position adjustment from the CIGAR string
	- Chromosome: The same for all three reads 
	- Position: The same for all three reads.  
	- UMI: The same for all three reads (from STL96.txt).  
	- Strandedness: The same for all three reads. 
	- Written: Only the first read should be written (Read 1) 
	
```
01:NS500451:154:HWKTMBGXX:1:11101:24260:1121:ACAGGACA	0	1	74	36	71M	*	0	0	TCCACCACAATCTTACCATCCTTCCTCCAGACCACATCGCGTTCTTTGTTCAACTCACAGCTCAAGTACAA	6AEEEEEEAEEAEEEEAAEEEEEEEEEAEEAEEAAEE<EEEEEEEEEAEEEEEEEAAEEAAAEAEEAEAE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
02:NS500451:154:HWKTMBGXX:1:11101:18996:1145:ACAGGACA	0	1	74	36	71M	*	0	0	GTCTCTTAGTTTATTATAAACCAGCTTCATAGGCCACAGAGGAAAAAGGACTATATACATACAGCCTTTTG	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:53G16	NH:i:1	HI:i:1	NM:i:2	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
03:NS500451:154:HWKTMBGXX:1:11101:25533:1187:ACAGGACA	0	1	74	36	71M	*	0	0	CTTGGTAACTTTCAGAGAATTAGTCACAACTTCTGAAGCAACCACAGTCCATGCAAGTCGACTGGTTTCTC	6AEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEAEEEEEEE<EEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```

4-6:
3 PCR duplicates which need position adjustment from the CIGAR string
	- Chromosome: The same for all three reads. 
	- Position: The same ADJUSTED pos for all three reads. The pos should be different, but corrections using CIGAR string should be the same.  
	- UMI: The same UMI from STL96.txt for all three reads. 
	- Strandedness: The same for all three reads. 
	- Written: Only the first read (04) should be written. 

```
04:NS500451:154:HWKTMBGXX:1:11101:6251:1098:AGTGCTGT	2	1	767	36	71M	*	0	0	GGCGTTCCAAACCACGGTCATCTCTTCTTTGCTTACTTTAGTGACTTCTGGAGGATCAGGGCGGCCAGGTC	/<EEAEEEEEEEEAEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
05:NS500451:154:HWKTMBGXX:1:11101:11995:1145:AGTGCTGT	2	1	768	36	1S70M	*	0	0	TGCATAACTCGTGCTGGTTTCCTCCTTTGTGGGGACGTGATAGGTCGAGTACCTGAAGTCTCTTCTTCTGT	6<EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
06:NS500451:154:HWKTMBGXX:1:11101:20566:1080:AGTGCTGT	2	1	769	36	2S69M	*	0	0	CCCCAAACAAGTGGTCCCCAAAAGAGACAGCCTCAAAATGGTCTAAGAAGCTGGCATAAAGGTCAGGAAAA	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEAEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```

7:
Unique read
	- Chromosome: Does not matter
	- Position: unique ADJUSTED position from all other reads in file 
	- UMI: any UMI from list STL96.txt
	- Strandedness: Does not matter
	- Written: Yes
	
```
07:NS500451:154:HWKTMBGXX:1:11101:23608:1152:GATCCTAG	0	1	77301	36	71M	*	0	0	GCCTTTCCGTTCCACCAGATAGTCTGTGATGATGCTCCCACCATCAAAGTCAGGCTTGGTCCAGCTCAGGA	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEE<AEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```

8-9:
2 Biological duplicates
	- Chromosome: Identical between the two 
	- Position: Identical between the two 
	- UMI: Different between the two, with both coming from the list STL96.txt
	- Strandedness: Same strandedness between the two reads.
	- Written: Yes, both are written.  

```
08:NS500451:154:HWKTMBGXX:1:11101:2638:1088:AGCTCTAG	0	1	767292	36	71M	*	0	0	CTCTGAGGTGGTTCCAGTAACATTCTTTAGTTCAAGCGTATATTCTCCCGTATCTCTTATAGTAGTTTCAC	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
09:NS500451:154:HWKTMBGXX:1:11101:1524:1089:AGTGCTGT	0	1	767292	36	71M	*	0	0	CTTTGCCAAAGAATGCTGTTTCTTTCTTTCTTTTCAACATGGAATCCAGTAACTTCTGAACCACCATCATA	6A/EE<AEEEAEEEEEEAEAEAE/AAEEEEEEEE/EEEEEAEEEEA6E<AEEEEAAE//EA66<A6EA<EA	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```
	
10-11:
2 Reads mapped to the same location (chrom, position, UMI identical) but the strandedness is different. 
	- Chromosome: Identical between the two 
	- Position: Identical between the two   
	- UMI: Identical between the two, chosen from the list STL96.txt
	- Strandedness: Different 
	- Written: Yes, both are written 

```
10:NS500451:154:HWKTMBGXX:1:11101:10568:1142:TCGTCATC	40	1	76750573	36	71M	*	0	0	TGCTGTCAGGATGTAAGCCCCACTATCCGTTCTTGAGGAATCCTTGTTTATCAGATGAGTGGAGAAGTCTG	6AEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
11:NS500451:154:HWKTMBGXX:1:11101:8096:1174:TCGTCATC	18	1	76750503	36	71M	*	0	0	CACCTGGCTTGATTTTCTGACCTGATTTATACCACGTCACACGAGGTTCTGGGTGGACGGTAATAGTTACT	/6<EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```
	
12-13:
2 Reads mapped to the same location (chrom, position) but the strandedness and UMIs are both different. 
	- Chromosome: Identical between the two 
	- Position: Identical between the two 
	- UMI: Different between the two, both from the list STL96.txt
	- Strandedness: Different between the two 
	- Written: Yes, both are written 

```
12:NS500451:154:HWKTMBGXX:1:11101:21621:1145:GGCGTATT	18	2	930	36	71M	*	0	0	TTCCACTGTTGCTTCATAACTGCAGTCCTAACATAAATGTCTGACATGTAGGATGATCTTAAGCAACCCCT	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEAAAEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
13:NS500451:154:HWKTMBGXX:1:11101:4845:1190:TCGACTTC	40	2	930	36	71M	*	0	0	CTTTGCCGACAGGGTTTTGCACTTCAAAGCTGTACACCCCACTGTCGCCAGGCACCACATTGATGATCTTG	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EAEEEEEEEEEAEEEEEAEEEEEEEAEEEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```

14-15:
2 Reads mapped to the same adjusted location (chrom, position, UMI identical) but the strandedness is different. 
	- Chromosome: Identical between the two reads.  
	- Position: Different position but identical ADJUSTED position between the two (use CIGAR string to specify soft-clipping for one read)  
	- UMI: Identical between the two (chosen from list STL96.txt)  
	- Strandedness: Different between the two: different Truth of (bflag & 16 ==16)  
	- Written: Yes, both are written.  

```
14:NS500451:154:HWKTMBGXX:1:11101:20178:1194:GCCGATAT	40	2	7671	36	5S66M	*	0	0	GCCTTTGATAGTGACAAATAGGCGCAAAGTGGCACTTGCACGCAGAGTGACCACCTTTCTGAGATCAGCAT	6AAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAEAEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
15:NS500451:154:HWKTMBGXX:1:11101:21370:1193:GCCGATAT	18	2	7595	36	71M	*	0	0	ATGACAGATACGGACTTGGTCCCACTAGCGTTCTTCACTGTCAGTGTGTATTTCCCACTGTCACTCCGGTC	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```

16:
Read mapped with incorrect UMI from given list: STL96.txt
	- Chromosome: Does not matter
	- Position: Unique from all other reads 
	- UMI: UMI which does not occur in the file STL96.txt 
	- Strandedness: Does not matter 
	- Written: Not if UMI list is specified in options, otherwise yes.  

```
16:NS500451:154:HWKTMBGXX:1:11101:25071:87858:TAGCATGG	0	2	7695	36	71M	*	0	0	GTTCCTCTTGCATTAATGGCTGTGACCCGATGTTCATATTCTAAGCCTTCAATAAGGCCGGTGGAGCGGTA	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```

17-18:
Two biological duplicates which differ from Pos but are corrected with Cigar string
	- Chromosome: The same for two reads.
	- Position: The pos should be different but ADJUSTED pos should be the same between the two. 
	- UMI: Two separate UMIs from STL96.txt
	- Strandedness: The same strand for both reads.  
	- Written: Yes, both are written. 

```
17:NS500451:154:HWKTMBGXX:1:11101:94095:71756:CGTTGGAT	40	2	76875	36	13S58M	*	0	0	GTGGGATGAGGCGCTCTTTTATATTGAGTTGGGCTGTGCAGGAGTCTTTTCCCACTTCATTGACGGCGTAG	6<EEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
18:NS500451:154:HWKTMBGXX:1:11101:69992:67325:CTAGAGGA	40	2	76879	36	17S54M	*	0	0	GTGGGATGAGGCGCTCTTTTATATTGAGTTGGGCTGTGCAGGAGTCTTTTCCCACTTCATTGACGGCGTAG	6<EEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
```

19-20:
Two reads which map to the exact same position but different strands.
	- Chromosome: The same for two reads.
	- Position: The pos should be different but ADJUSTED pos should be the same between the two. 
	- UMI: Two separate UMIs from STL96.txt
	- Strandedness: Reads come from different strands. 
	- Written: Yes, both are written. 

```
19:NS500451:154:HWKTMBGXX:1:11101:10266:1114:TTCGTTCG	40	3	768	36	4S67M	*	0	0	GTGAAACTCGGTGGGATGAGGCGCTCTTTTATATTGAGTTGGGCTGTGCAGGAGTCTTTTCCCACTTCATT	6<EAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEAEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
20:NS500451:154:HWKTMBGXX:1:11101:4191:1194:TTCGCCTA	18	3	693	36	10M10D40M11S	*	0	0	GTTCCTCTTGCATTAATGGCTGTGACCCGATGTTCATATTCTAAGCCTTCAATAAGGCCGGTGGAGCGGTA	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```

21-24:
Four reads mapped to the same position, two are biological duplicates with one PCR duplicate per.  
	- Chromosome: The same for all four reads.
	- Position: The pos can be different but ADJUSTED pos should be the same between the two. 
	- UMI: Two separate UMIs from STL96.txt, twice replicated each.
	- Strandedness: Reads come from the same strand. 
	- Written: Yes, both biological replicates are written, however the PCR duplicates are not. (Read 21 and 22 are written) 

```
21:NS500451:154:HWKTMBGXX:1:11101:53424:50433:CAAGTCGT	0	3	7877	36	2S69M	*	0	0	GTGAAACTCGGTGGGATGAGGCGCTCTTTTATATTGAGTTGGGCTGTGCAGGAGTCTTTTCCCACTTCATT	6<EAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEAEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
22:NS500451:154:HWKTMBGXX:1:11101:347:19169:TATAGCGC	0	3	7890	36	15S56M	*	0	0	GTGGGATGAGGCGCTCTTTTATATTGAGTTGGGCTGTGCAGGAGTCTTTTCCCACTTCATTGACGGCGTAG	6<EEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
23:NS500451:154:HWKTMBGXX:1:11101:6610:1112:CAAGTCGT	0	3	7880	36	5S66M	*	0	0	GTGGGATGAGGCGCTCTTTTATATTGAGTTGGGCTGTGCAGGAGTCTTTTCCCACTTCATTGACGGCGTAG	6<EEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
24:NS500451:154:HWKTMBGXX:1:11101:22058:1205:TATAGCGC	0	3	7875	36	71M	*	0	0	CTTTGGTCCAGGCAACCTCTGGGAATGGTTTGCCTCTCACACCTGCCTCGAGTCTAATGGTGTCCCCTGCT	6AEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```

25-28:
This is a tricky one. It is actualy two separate mappings. 4 reads: I want to split up the lines of a mapping of biological replicate. R1 is kept, R2 is a whole new mapping, R3 is biological duplicate, R4 is a PCR duplicate. 
	- Chromosome: The same for all four reads.
	- Position: The pos MUST be different but ADJUSTED pos should be the same between the first one  and the last two. 
	- UMI: Different UMIs all of which are from the list STL96.txt.
	- Strandedness: Reads come from the same strand. 
	- Written: Yes, both biological replicates are written, as well as the novel sequence between them however the PCR duplicate is not (Reads 25, 26, 27)

```
25:NS500451:154:HWKTMBGXX:1:11101:23713:1207:AAGGTACG	0	3	76200105	36	5S66M	*	0	0	GGCTTTGCAGGTGTACTCCCCGCAGTCAACGACCTGGGTTCTCAGGATCTCCAGGCTGGAGACATAGTTCG	/AEEEEEEEEEEEE/EEEEEEEEEEEE/EEEAEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEAEEA<EEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
26:NS500451:154:HWKTMBGXX:1:11101:6943:1217:AGCTACCA	0	3	76200105	36	6S65M	*	0	0	CACTTTGCTAGGCTTGCCAATGCCCACAATGTTCTCGGCAGAAACCCTGAACTCATATTCAAGCCCTTCAT	/AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
27:NS500451:154:HWKTMBGXX:1:11101:5571:1221:AGGACATG	0	3	76200110	36	10S61M	*	0	0	TTCCAGGTACACAAAAGTCTTCTGAGTAAACAACCTGTACTTTTTGCTACTTCGGATCTGCTTCTTGTCTT	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
28:NS500451:154:HWKTMBGXX:1:11101:86637:67646:AGGACATG	0	3	76200120	36	20S51M	*	0	0	TTCCAGGTACACAAAAGTCTTCTGAGTAAACAACCTGTACTTTTTGCTACTTCGGATCTGCTTCTTGTCTT	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
```















Dont need the rest of these reads I think

```
NS500451:154:HWKTMBGXX:1:11101:89435:58226:TCGTAGGT	0	2	76708822	36	71M	*	0	0	CCTGCTCATGGATACTCTTGAAGGATTGTCCAGTAAACTGGAAGCTGGTTTTGGAAGTTCCTCCTTCACCA	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:20948:1192:ACTGTCAG	0	2	76764786	36	71M	*	0	0	CACCTTTTCCAGCCTTCGTCAGGTGATGCGTCTGCTATTTTGGGTCTCATCTCTACGACATATCCAATTAT	6<EEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEAEEEEAEEEEEEA	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:19377:1220:TCGTAGGT	0	2	76708822	36	71M	*	0	0	CCTGCTCATGGATACTCTTGAAGGATTGTCCAGTAAACTGGAAGCTGGTTTTGGAAGTTCCTCCTTCACCA	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:4794:96445:ACTGTCAG	0	2	76764786	36	71M	*	0	0	CACCTTTTCCAGCCTTCGTCAGGTGATGCGTCTGCTATTTTGGGTCTCATCTCTACGACATATCCAATTAT	6<EEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<EEAEEEEAEEEEEEA	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:13631:1222:TCGTAGGT	0	2	76944432	36	71M	*	0	0	GCCAGTACCAGTAGATTCGGTCAGATCGCTCAATTTTTACACCATTTTTATACCATTCACACTCAGGGTCT	/AEEEEEEEEEEEEEEEEEEEEEEAEEE<EEEEEEEEEEE6<EEEEEEEEEEEEEEEEEEAEEEEE<AEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:12354:25636:TCGTAGGT	0	2	76708822	36	71M	*	0	0	CCTGCTCATGGATACTCTTGAAGGATTGTCCAGTAAACTGGAAGCTGGTTTTGGAAGTTCCTCCTTCACCA	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:5653:1249:CAACGTTG	0	2	52287124	36	71M	*	0	0	GTTTGTAGTTGACGTTGGTAGCCACGTCCTGGGCCATCCTTGCAGCAGTGATGCTGACCATGTCTCCAGGG	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEAEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:38067:55266:TAGCAAGG	0	2	76718924	36	71M	*	0	0	GTTCCTCTTGCATTAATGGCTGTGACCCGATGTTCATATTCTAAGCCTTCAATAAGGCCGGTGGAGCGGTA	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:13611:1237:GACATGAG	0	2	76748739	36	71M	*	0	0	GTGGGATAACAAGGTCTTTAGCAATTACTGGCACACTCAGTTGTCTGGGGTCACTGATGCCCTTTTCATTC	6AEEEEEEEEEEEEEEEEEEAAEA/EEAEEAE<EEEE/EE/EEEEEEE<E/AA<EA6E/6A6AAEEEE</E	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:96518:49641:CACACACA	0	2	76723334	36	71M	*	0	0	CCACGATCACCCAATTGAGCCTGCTAGTCTCCCGTCTCTCTACAATATAGTGAGTTATCTCTGCTCCTCCA	6/EEEEEEAAEEEEEE<EAEEAEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEAEEEEEEEEEEEEEEEAEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:12146:1234:CACACACA	0	2	76723334	36	71M	*	0	0	CCACGATCACCCAATTGAGCCTGCTAGTCTCCCGTCTCTCTACAATATAGTGAGTTATCTCTGCTCCTCCA	6/EEEEEEAAEEEEEE<EAEEAEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEAEEEEEEEEEEEEEEEAEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:6311:1206:TCGTCATC	0	2	76742000	36	71M	*	0	0	GTTGCTGTGCCTGCACTGTTTGTTGCTGTGACACTGTATTTCCCGAAGTCATCTTTGCTGCTTTCTTTAAT	6<EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEAEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:14374:1195:TGAGTGAG	0	2	76744713	36	71M	*	0	0	CTTCTTCAAGGTTCAGTGCTTTAAACTGGGTGTCATGAATAATGGTTTTGTTGACTTTTGTCCACAAAATA	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:13546:46364:TCAGGACT	0	2	76737558	36	71M	*	0	0	CTGCCATAGTCCGAGTTCCTGCCACATTCTTCAGTGTAAGGACATACTGCCCAGTGTCTCTTCTTACACAG	6<AAEEEEEEEEEA<A/EEEEEAEEEEEEEEEEEEEEEEEE/AEE/EEE<EAEAEEEEAA/EEEEEEEEAE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:18102:1218:TCAGGACT	0	2	76737558	36	71M	*	0	0	CTGCCATAGTCCGAGTTCCTGCCACATTCTTCAGTGTAAGGACATACTGCCCAGTGTCTCTTCTTACACAG	6<AAEEEEEEEEEA<A/EEEEEAEEEEEEEEEEEEEEEEEE/AEE/EEE<EAEAEEEEAA/EEEEEEEEAE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:9862:1266:GCCGATTA	0	2	76736951	36	71M	*	0	0	GCTCGAACCAAGGGAGAAGTTTCACTTGGAAGACTCAGGCCGGCAGCATTTTCAGCAAAAACCCGATACTC	6AEEEA/EEEEAA/EEEEAEAEEEEE/EEEEEEAEEAAAEEEAEEEEAEAEEEEE<EEE/EEE6EEEEE/<	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:17732:1263:CGCCTATT	0	2	30402279	36	71M	*	0	0	GGCTCCTGCTGTGAGAAGGGATGGGCTGGACTGGCTCTCTTCTCCATGGTTCATGGAGTGGATGATACATT	6<EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:22329:8818:TCAGGACT	0	2	76737558	36	71M	*	0	0	CTGCCATAGTCCGAGTTCCTGCCACATTCTTCAGTGTAAGGACATACTGCCCAGTGTCTCTTCTTACACAG	6<AAEEEEEEEEEA<A/EEEEEAEEEEEEEEEEEEEEEEEE/AEE/EEE<EAEAEEEEAA/EEEEEEEEAE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:15723:1090:AGACACTC	0	2	52169122	36	42M679N29M	*	0	0	GTAGTGGCCTTTGCTCTTCTCATAGTTCTTCTTGTACTGGTAAGAACTGGCGAGCTGGCTTGTATTTAGGA	/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEA	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
NS500451:154:HWKTMBGXX:1:11101:1232:1273:GGATAACG	0	2	52308305	36	71M	*	0	0	CTTCGTATTCCTGGGTAATGGTCTGGGGGAAGAAGCCTTTGCCTCTGTCCTCTTCATACTCAGCTTTGTAG	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<<EEEEEEEEEEA<EEEEEEEEE<EEEEAAEE<EEEAEA	MD:Z:4A66	NH:i:1	HI:i:1	NM:i:1	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:21274:1224:CTAGGAAG	0	2	76771715	36	71M	*	0	0	CCCGCCATCATACTTGGGAGGGTTCCAAGTCAAAGACACAGTATTTTTAGTGACTTCTTTGACTTCCAGGT	6AEEEEEEEEEEEEEEEEEEEEEEAEAAAEEEEEEEEAEEEEEEEEEEEAEEEEAEEEEEE<AEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:26383:1289:CAACGATC	0	2	76797980	36	71M	*	0	0	GTTGGGGCTTCTATGATGAGCAGCTCTGCCACGGATTTATCTTGCCCAGCAGTGACGACGTATTCACCTTC	6//E/AEEAEE6EA<EEEEEEEEE/E//EAEE//EEEAEEE</<AE<EEEE/EE<AE<A/EE/E<AEEE/6	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:9957:51931:TAGCAAGG	0	2	76718924	36	71M	*	0	0	GTTCCTCTTGCATTAATGGCTGTGACCCGATGTTCATATTCTAAGCCTTCAATAAGGCCGGTGGAGCGGTA	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:8042:1317:GCGGTATT	0	2	76747236	36	71M	*	0	0	CCCCTGATATAGCAACAGGTCCCTCAGGTGGCCCTGGTCTATCGAGAACTTTGACATTGATAGTAACTGAT	/A/EEEEEE6EEAEEEEAEEEEAEEAAAE<E<EAEEAEEEAEE<E6EEEEEEAAAA/E<AEEE</AA/A6E	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:52657:6101:TAGCAAGG	0	2	76718924	36	71M	*	0	0	GTTCCTCTTGCATTAATGGCTGTGACCCGATGTTCATATTCTAAGCCTTCAATAAGGCCGGTGGAGCGGTA	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:59531:27832:TCAGGACT	0	2	76737558	36	71M	*	0	0	CTGCCATAGTCCGAGTTCCTGCCACATTCTTCAGTGTAAGGACATACTGCCCAGTGTCTCTTCTTACACAG	6<AAEEEEEEEEEA<A/EEEEEAEEEEEEEEEEEEEEEEEE/AEE/EEE<EAEAEEEEAA/EEEEEEEEAE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:14812:1379:CGCCTATT	0	2	76750363	36	71M	*	0	0	GGCAACATTTGTCCATGCCAACCTGCTGGTTTCACGTTTCTGTACGATGTAATGATCAATTTTGGCACCTC	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:14449:1105:GAACAGGT	0	2	76821131	36	54M1266N17M	*	0	0	GCTCCACTTTTTCCACAACCTCGGGTTTTTCAGGAACTTTCTTCTTTGAAACAGCTTTAAGTTGGAACATT	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
NS500451:154:HWKTMBGXX:1:11101:5808:1335:TCGACTTC	0	2	76798066	36	71M	*	0	0	GTCCTTGATGATTAGAGAGTGTTTGTACTTATCAATTCGGTATGATATACGGTTGTCAAAAGCCACTTCTT	6AEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEAEAEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:22945:1315:AGTGCTGT	0	2	76776807	36	71M	*	0	0	GGTGTCATAAACAACAGGCTCTTTGCTGTCGGGTTTCTTTGGCGGAGCCTTGAACCAGGTGAGAGTTGGGA	66EEEEEEEE6EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:5447:1346:ATCGTTGG	0	2	76770826	36	71M	*	0	0	CTTCAAGCAGTTTAGGGACCTTGCATGTAGTTTTGGCACTGGCAGATGTCACAGGCATCCAGATGTCTTTA	6EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:2312:1365:TCCACATC	0	2	76942250	36	71M	*	0	0	GTTCGGTTAATCACAAGTCGCTGTTTGGTGCCTTTGACAATGGATTGCACACGATCATCAGGCTTGATTTG	6AEEEEEEEEEEEEEEEEEE<EEEEEEEEEEEEEEAEEEEEAEEEEEEEEEEEEEEEEEEEAEEEEEEEE/	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:22228:1353:CACACACA	0	2	76745256	36	71M	*	0	0	CCACCACATTGGAAACTGTGATTCCATACTGTCCACCATCATCCTTATGAGTTTCTTTAATGCTGAGTGTG	6<EEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:22955:1351:ATCGAACC	0	2	76887211	36	71M	*	0	0	CTCGTGTCAAAATTTTTGGAGGTGATTTTAAATTTCTTGCTGCTTCGCACTTGCTTCCGATCTTTAACCCA	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEA6A6	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:13771:2468:GGAATTGG	0	2	60763291	36	52M16388N19M	*	0	0	CTCTTTTTCTGTCCTCCATCCGCAAACTTGCACAGTAAAGGTTCTGTAGGAGCAGAAACTCCTGGTGGGGT	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
NS500451:154:HWKTMBGXX:1:11101:63966:22359:AGCTACCA	0	2	76814762	36	71M	*	0	0	TGTCAAAACGATGTCTTTGATTTCTTTCACGAATTTCAGTGGCACGGCCTTGAGTTGGTAAGTGAACGGAG	66EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:42688:27141:GGAATTGG	0	2	60763291	36	52M16388N19M	*	0	0	CTCTTTTTCTGTCCTCCATCCGCAAACTTGCACAGTAAAGGTTCTGTAGGAGCAGAAACTCCTGGTGGGGT	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
NS500451:154:HWKTMBGXX:1:11101:5434:1427:GGAATTGG	0	2	60763291	36	52M16388N19M	*	0	0	CTCTTTTTCTGTCCTCCATCCGCAAACTTGCACAGTAAAGGTTCTGTAGGAGCAGAAACTCCTGGTGGGGT	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
NS500451:154:HWKTMBGXX:1:11101:17250:1435:AGCTACCA	0	2	76814762	36	71M	*	0	0	TGTCAAAACGATGTCTTTGATTTCTTTCACGAATTTCAGTGGCACGGCCTTGAGTTGGTAAGTGAACGGAG	66EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:52106:17419:AGCTACCA	0	2	76814762	36	71M	*	0	0	TGTCAAAACGATGTCTTTGATTTCTTTCACGAATTTCAGTGGCACGGCCTTGAGTTGGTAAGTGAACGGAG	66EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:1646:11667:AACGCCAT	0	2	76749268	36	71M	*	0	0	CGTCAATGTGACAGATGTCTTTGTGACCTCCTTTATCTTTAGATTCTGTGGAGGGCCTGGCGTATCTAGCA	/AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA<	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:26287:1436:AACGCCAT	0	2	76749268	36	71M	*	0	0	CGTCAATGTGACAGATGTCTTTGTGACCTCCTTTATCTTTAGATTCTGTGGAGGGCCTGGCGTATCTAGCA	/AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA<	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:5518:1226:GGATAACG	0	2	52328864	36	55M1890N16M	*	0	0	ATTGGCCATAGTCCGATGTCCTTGTCGTTGTGGTCTCATAGACTTCTGTAATTGTCTCTCCTGGCACCTCT	6EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEAEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
NS500451:154:HWKTMBGXX:1:11101:2260:1448:CTAGGAAG	0	2	76902958	36	71M	*	0	0	GGGGTGTGGATGGTACAGCAACCAGAGAATGTTCACAGAAGAGAGAGACCCTTCTAAGGCTTCCTTTCAGA	/<AE/AEEEEEEAEAE<EEAEEEEAEEEEAEEEEEEEEEEEEEEEEE6EEEEEAEEAEEEAEEEEEEEAEA	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:24998:1299:GTTCTGCT	0	2	76784182	36	36M103N35M	*	0	0	CAACCGGTTTGAGGTCAAGAACTGGGCCAGGGACATCGAACACTTCAACCCGGACTGCGGCAAATTTGGAA	66EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
NS500451:154:HWKTMBGXX:1:11101:6968:1257:GAGAAGTC	0	2	76788807	36	13M83N58M	*	0	0	GGTGAAGGAGGACTTGTTGGGTCTTCGATTGACAGAATTTCTGTGGGCTCACTTGGGTGACCAACTCCTGC	6AEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEE<EAEEEEEEA	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
NS500451:154:HWKTMBGXX:1:11101:14739:1447:ACGTCAAC	0	2	76730808	36	71M	*	0	0	GGAATCAGTAGTTGGACATCTTCACCAGCTTTAGCTATAACAGATGTTCTGAGAGCCACATCTAGGTCGAT	//EEEEE6EEA/AEAEEEE<AEAE/AE<EAEEEE/EEEEEEEEE//EEEAEEEEA/AA/EE/AEA/EEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
NS500451:154:HWKTMBGXX:1:11101:24936:1293:TGAGTGAG	0	2	52159545	36	25M1086N46M	*	0	0	GTCCGATTGCTTCTTTATACTTGACCGAGCTAATGTGGTCCTGCGTTTGCTTCACTCTGAGCATCTCAGGC	6AEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU	XS:A:-	XG:Z:A
```
