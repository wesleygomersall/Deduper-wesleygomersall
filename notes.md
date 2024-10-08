
# be careful wiht CIGAR String from sam file. Due to soft clipping the reads may have a different position in the alignment however they are in fact PCR duplicates.

# After looking at the positions in the alignment output, check UMI barcodes. PCR duplicates should contain the same UMI sequence. 

# If two reads are PCR duplicates: 
# they will share : 
#        Chromosome
#        5' start position of read
#        strand (+ or -) 
#        UMI 

# What is the difference between single-end vs paired-end.

# This removal tool will be for single-end data only. 

# ~~~~~~~~~~~~~~~~~

# Given a SAM file of Uniquely mapped reads, remobe all PCR duplicates (retain only a dinlge copy of each read) 
# do not need to sort
# do not need to check if reads are mapped or not (they shold all be mapped) 
# need to adjust for soft clipping
# single-end reads
# UMIs are known
# considerations: 
#     millions of reads - aboid loading everything into memory!
#     be sure to utilize functions appropriately 
# CHALLENGE: 
#     single-end vs paired end
