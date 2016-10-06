# -*- coding: utf-8 -*-
from scipy import stats
import numpy as np
from sys import argv,exit
import argparse
import pysam
from Bio import SeqIO
import pandas
from collections import defaultdict
from math import log

script,genbankfile,tss = argv

def parse_genbank(genbankfile):
	try:
		genome =SeqIO.read(open(genbankfile,"rU"),"genbank")
	except ValueError: #enthält die Datei mehrere Genomsequenzen wird das Programm beendet.
		exit("RefSeq-File contains multiple sequences. Please provide whole genome sequence!!!\nProgram terminated!")
	array_of_features_plus = np.zeros(len(genome.seq),dtype=int) #Array of Dicts für die Features
	array_of_features_minus = np.zeros(len(genome.seq),dtype=int)
	for feature in SeqIO.read(open(genbankfile,"rU"),"genbank").features:
		if feature.type == "source": #"Source", also die gesamte Referenzsequenz wird als ein normales Feature behandelt, was für Probleme sorgt
			next
		elif feature.location.strand == 1:
			for i in range(int(feature.location.start),int(feature.location.end)): 
				array_of_features_plus[i] = 1
		elif feature.location.strand == -1:
			for i in range(int(feature.location.start),int(feature.location.end)):
				try:
					array_of_features_minus[i] = 1
				except IndexError:
					pass
	return array_of_features_plus,array_of_features_minus


def parse_csv(tss):
	df = pandas.read_csv(tss,sep=",")
	df.columns = [c.replace(" ","_") for c in df.columns]
	tss = df.Position.tolist()
	strang = df.Strand.tolist()
	reads = df.No_Read_Starts.tolist()
	return tss,strang,reads

def classify_tss(plus,minus,tss,strang,reads):
	for i,j,k in zip(tss,strang,reads):
		if j == "Fwd":
			if plus[i] == 1:
				print i,j,k,"iTSS"
			elif plus[i] == 0:
				for n in plus[i:i+301]:
					if n == 1:
						print i,j,k,"pTSS"
						break
					else:
						pass
				else:
					for n in  minus[i-100:i+100]:
						if n == 1:
							print i,j,k,"asTSS"
							break
					else:
						print i,j,k,"oTSS"
		elif j =="Rev":
			if minus[i] == 1:
				print i,j,k,"iTSS"
			elif minus[i] == 0:
				for n in minus[i-301:i]:
					if n == 1:
						print i,j,k,"pTSS"
						break
					else:
						pass
				else:
					for n in  plus[i-100:i+100]:
						if n ==1:
							print i,j,k,"asTSS"
							break
					else:
						print i,j,k,"oTSS"


plus,minus = parse_genbank(genbankfile)
tss,strang,reads = parse_csv(tss)
classify_tss(plus,minus,tss,strang,reads)
