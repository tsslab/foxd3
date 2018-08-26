#!/usr/bin/env python
from Bio import SeqIO
import sys

read_length = int(sys.argv[2])
quality=read_length*'H'

out = open(sys.argv[1].split('.')[0]+"_"+sys.argv[2]+"bp_reads.fastq", "w")

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
	chrom = str(seq_record.id)
	seq = str(seq_record.seq).upper()
	print chrom, len(seq)
	#print seq
	for i in range(len(seq) - read_length + 1):
		#print i
		out.write("@synthetic_read-"+sys.argv[1].split('.')[0]+':'+chrom+':'+str(i)+':'+str(i+read_length)+'\n')
		read=str(seq[i:i+read_length])
		#print len(read)
		#print read
		out.write(read+'\n')
		out.write("+\n")
		out.write(quality+'\n')
	
