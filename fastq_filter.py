# -*- coding: utf-8 -*-
from Bio import SeqIO
import re
from sys import argv

script,fastq = argv

fastq_new = str(fastq)+"_filtered.fastq"
outfile = open(fastq_new,"w")
counter = 0
for record in SeqIO.parse(fastq,"fastq"):
	a_counter = 0
	for char in str(record.seq):
		if a_counter > 1:
			break
		elif char == "A":
			a_counter += 1
	else:
		SeqIO.write(record,outfile,"fastq")
		counter += 1
outfile.close()
print "%d Reads enthalten kein A" %counter
