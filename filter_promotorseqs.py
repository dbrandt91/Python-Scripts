# -*- coding: utf-8 -*-
import numpy as np
from sys import argv,exit
import pysam
from Bio import SeqIO
import pandas

script,tss = argv

with open(tss,"r") as f:
	for count,seq in enumerate(f):
		for l_count,letter in enumerate(seq.rstrip()):
			if letter.isupper():
				if seq[l_count:l_count+2] != "TG":
					print ">%d\n%s" %(count,seq.rstrip())
				break
