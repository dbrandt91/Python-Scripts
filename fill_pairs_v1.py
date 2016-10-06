# -*- coding: utf-8 -*-
from sys import argv,exit
import argparse
import pysam
from Bio import SeqIO
import subprocess

script, bamfile, fasta = argv

index_of_dot = bamfile.index(".bam")
filename_out = bamfile[:index_of_dot]+"_filled_pairs.fasta"
bamfile_sorted = bamfile[:index_of_dot]+"_namesorted.bam"
pysam.sort("-n",bamfile,"-@",3,"-o",bamfile_sorted)
outfile = open(filename_out, "w")
fasta = open(fasta,"rU")
start = 0
end = 0
name = ""
i = 0
read_tuples = []
f = pysam.AlignmentFile(bamfile_sorted,"rb")
for read in f:
	if read.is_proper_pair and read.is_read1 and not read.is_reverse and not read.is_duplicate:
		mate = f.next()
		if not mate.is_unmapped:
			read_tuples.append((str(read.query_name),int(read.reference_start),int(mate.reference_end)))
	elif read.is_proper_pair and read.is_read1 and read.is_reverse and not read.is_duplicate:
		mate = f.next()
		if not mate.is_unmapped:
			read_tuples.append((str(read.query_name),int(read.reference_end),int(mate.reference_start)))
	elif not read.is_proper_pair and not read.is_reverse:
		read_tuples.append((str(read.query_name),int(read.reference_start),int(read.reference_end)))
	elif not read.is_proper_pair and read.is_reverse:
		read_tuples.append((str(read.query_name),int(read.reference_end),int(read.reference_start)))
	elif read.is_unmapped:
		pass
for record in SeqIO.parse(fasta,"fasta"):
    for read in read_tuples:
        if read[1] < read[2]:
            #outfile.write(">%s\n%s\n" %(read[0],str(record.seq[int(read[1])-1:int(read[2])]).upper()))
            print "%s,%d,%d" %(read[0],int(read[1])-1,int(read[2]))
        elif read[1] > read[2]:
            #outfile.write(">%s\n%s\n" %(read[0],str(record.seq[int(read[2])-1:int(read[1])].reverse_complement()).upper()))
            print "%s,%d,%d" %(read[0],int(read[2])-1,int(read[1]))

outfile.close()
fasta.close()
