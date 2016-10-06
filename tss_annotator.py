from Bio import SeqIO
import pandas
from collections import defaultdict
import numpy as np
from sys import exit,argv

script,genbankfile,tss = argv

def read_tssfile(tss):
	df = pandas.read_csv(tss,sep=",")
	tss_array = df.TSS.tolist()
	strand_array = df.Strand.tolist()
	return tss_array,strand_array

def parse_genbank(genbankfile):
	try:
		genome =SeqIO.read(open(genbankfile,"rU"),"genbank")
	except ValueError:
		exit("RefSeq-File contains multiple sequences. Please provide whole genome sequence!!!\nProgram terminated!")
	array_of_features_plus = np.zeros(len(genome.seq),dtype=dict)
	array_of_features_minus = np.zeros(len(genome.seq),dtype=dict)
	for feature in SeqIO.read(open(genbankfile,"rU"),"genbank").features:
		if feature.type == "source":
			next
		elif feature.location.strand == 1:
			for i in range(int(feature.location.start)-300,int(feature.location.end)):
				array_of_features_plus[i] = feature.qualifiers.copy()
				array_of_features_plus[i].update({"strand":"Fwd"})
		elif feature.location.strand == -1:
			for i in range(int(feature.location.start),int(feature.location.end)+300):
				try:
					array_of_features_minus[i] = feature.qualifiers.copy()
					array_of_features_minus[i].update({"strand":"Rev"})
				except IndexError:
					pass
	return array_of_features_plus,array_of_features_minus
	

def get_start(array_of_features_plus,array_of_features_minus,tss_array,strand_array):
	tss_start_dict = defaultdict(lambda:0)
	for tss,strand in zip(tss_array,strand_array):
		try:
			if  not array_of_features_plus[tss] and not array_of_features_minus[tss]:
				pass
			elif array_of_features_plus[tss] and strand == "Fwd":
				print "%d,%s,%s"%(tss,str(array_of_features_plus[tss]["locus_tag"]).translate(None,"['']") ,str(array_of_features_plus[tss]["gene"]).translate(None,"['']") )
			elif array_of_features_minus[tss] and strand =="Rev":
				print "%d,%s,%s" %(tss,str(array_of_features_minus[tss]["locus_tag"]).translate(None,"['']") ,str(array_of_features_minus[tss]["gene"]).translate(None,"['']") )
		except KeyError:
			pass

tss_array,strand_array = read_tssfile(tss)
array_of_features_plus,array_of_features_minus = parse_genbank(genbankfile)
get_start(array_of_features_plus,array_of_features_minus,tss_array,strand_array)
