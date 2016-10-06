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

#BAM-File einlesen und Coverage-Array berechnen
def read_and_parse_bam_to_array(bamfile,genome):
	genome_length = len(genome.seq)
	coverage_plus = np.ones(genome_length, dtype=int) #Read-Coverage wird in einem numpy-Array mit der Länge des Genoms berechnet
	coverage_minus = np.ones(genome_length, dtype=int) #strangspezifisch
	for read in pysam.AlignmentFile(bamfile,"rb"):
		if read.is_reverse == False:
			coverage_plus[read.reference_start:read.reference_start+read.query_length] += 1
		elif read.is_reverse == True:
			coverage_minus[read.reference_start:read.reference_start+read.query_length] += 1
		else:
			next
	return coverage_plus, coverage_minus

#Funktion zum Parsen von Genbank-Files (nicht schön aber selten)
def parse_genbank(genbankfile):
	try:
		genome =SeqIO.read(open(genbankfile,"rU"),"genbank")
	except ValueError: #enthält die Datei mehrere Genomsequenzen wird das Programm beendet.
		exit("RefSeq-File contains multiple sequences. Please provide whole genome sequence!!!\nProgram terminated!")
	array_of_features_plus = np.zeros(len(genome.seq),dtype=dict) #Array of Dicts für die Features
	array_of_features_minus = np.zeros(len(genome.seq),dtype=dict)
	for feature in SeqIO.read(open(genbankfile,"rU"),"genbank").features:
		if feature.type == "source": #"Source", also die gesamte Referenzsequenz wird als ein normales Feature behandelt, was für Probleme sorgt
			next
		elif feature.location.strand == 1:#Feature-Dicts auf dem Forward-Strang werden an alle Positionen, für die sie annotiert sind in den Array kopiert (+300 bp UTR)
			for i in range(int(feature.location.start)-300,int(feature.location.end)): 
				array_of_features_plus[i] = feature.qualifiers.copy()
				array_of_features_plus[i].update({"strand":"fwd"})
		elif feature.location.strand == -1:
			for i in range(int(feature.location.start),int(feature.location.end)+100): #dasselbe für den Reverse-Strang (nur 100 bp UTR)
				try:
					array_of_features_minus[i] = feature.qualifiers.copy()
					array_of_features_minus[i].update({"strand":"rev"})
				except IndexError:
					pass
	return array_of_features_plus,array_of_features_minus,genome
	
#TSS-Erkennung auf dem forward-Strang
def tss_detection_plus(coverage_array,window_size,step_height,min_p_value,min_slope):
	mean_coverage= sum(coverage_array)/len(coverage_array) #Berechnung der mittleren Coverage
	tss_list = [] #leere Liste für die detektierten TSS
	k = 1 #entspricht der Position im Genom
	counter = 0  #Counter um direkt angrenzende TSS (< 5 bp) auszuschließen
	for j,i in enumerate(np.diff(coverage_array)): #i entspricht der Position im Differenzvektor
		if counter > 0: #Ist counter größer als 0 wird er immer um 1 dekrementiert
			counter -= 1
			k+= 1 
			
		elif i >= step_height and float(coverage_array[j+1])/float(coverage_array[j]) > 1.3: #alle Stufen über 16 Reads werden direkt gezählt
			tss_list.append(k+1)
			k += 1
			counter = 5
			
		elif mean_coverage <= i < step_height: #Stufen zwischen der mittleren Coverage und 16 Reads werden gezählt, wenn die Regression signifikant ist
			slope,intercept,r_value,p_value,std_err=stats.linregress(np.arange(k-(window_size)/2,k+(window_size+2)/2),coverage_array[np.arange(k-(window_size)/2,k+(window_size+2)/2)])
			if p_value < min_p_value and slope > min_slope:
				tss_list.append(k+1)
				k += 1
				counter = 5
			else:
				k+=1
				
		else:
			k+=1
	tss_array = np.array(tss_list,dtype=int) #Umwandeln der TSS-Liste in numpy-Array
	print len(tss_array)
	return tss_array

#TSS-Erkennung (-)-Strang. Funktioniert wie auf dem fwd-Strang, nur wird der Coverage-Array rückwärts durchlaufen und die Indizes entsprechend invertiert.
def tss_detection_minus(coverage_array,window_size,step_height,min_p_value,min_slope):
	mean_coverage= sum(coverage_array)/len(coverage_array)
	tss_list = []
	k = len(coverage_array)
	counter = 0
	for j,i in enumerate(np.diff(coverage_array[::-1])):
		if counter > 0:
			counter -= 1
			k-= 1
		elif i >= step_height and float(coverage_array[::-1][j+1])/float(coverage_array[::1][j]) > 1.3:
			tss_list.append(k-1)
			k -= 1
			counter = 5
		elif mean_coverage <= i < step_height:
			slope,intercept,r_value,p_value,std_err=stats.linregress(np.arange(k-(window_size)/2,k+(window_size+2)/2),coverage_array[np.arange(k-(window_size)/2,k+(window_size+2)/2)])
			if p_value < min_p_value and slope > min_slope:
				tss_list.append(k-1)
				k -= 1
				counter = 5
			else:
				k-=1
		else:
			k-=1
	tss_array = np.array(tss_list,dtype=int)
	print len(tss_array)
	return tss_array

#Zuordnung der Position der annotierten TSS zu einem Feature (wenn möglich)
def get_start_plus(array_of_features,tss_array,genome):
	seq = str(genome.seq) #Genomsequenz
	tss_start_dict = defaultdict(lambda:0) #leeres Dict mit Standardwert 0
	for tss in tss_array:
		if  not array_of_features[tss]: #Wenn Feature noch nicht annotiert
			array_of_features[tss]= {"Novel Transcript":"Yes"} #Neuannotation des Transkriptes
			array_of_features[tss].update({"strand":"fwd"}) 
			array_of_features[tss].update({"Sequenz(-59-+3)":seq[tss-59:tss+2]}) #Genomsequenz im Bereich der TSS
			tss_start_dict[tss] = array_of_features[tss].copy()
		elif array_of_features[tss] and array_of_features[tss]["strand"] == "fwd": #Für bereits annotierte Features
			array_of_features[tss].pop("translation",None) #wird das Feature-Dict aus der Genbank-File etwas beschnitten...
			array_of_features[tss].pop("transl_table",None)
			array_of_features[tss].pop("organism",None)
			array_of_features[tss].update({"Novel Transcript":"No"})
			array_of_features[tss].update({"Sequenz(-59-+3)":seq[tss-59:tss+2]})
			tss_start_dict[tss] = array_of_features[tss].copy()#... und dann kopiert
			
	return tss_start_dict

def get_start_minus(array_of_features,tss_array,genome):
	comp_seq = str(genome.seq.complement())
	tss_start_dict = defaultdict(lambda:0)
	for tss in tss_array:
		if  not array_of_features[tss]:
			array_of_features[tss]={"Novel Transcript":"Yes"}
			array_of_features[tss].update({"strand":"rev"})
			array_of_features[tss].update({"Sequenz(-59-+3)":comp_seq[tss-3:tss+58][::-1]})
			tss_start_dict[tss] = array_of_features[tss].copy()
		elif array_of_features[tss] and array_of_features[tss]["strand"] == "rev":
			array_of_features[tss].pop("translation",None)
			array_of_features[tss].pop("transl_table",None)
			array_of_features[tss].pop("organism",None)
			array_of_features[tss].update({"Novel Transcript":"No"})
			array_of_features[tss].update({"Sequenz(-59-+3)":comp_seq[tss-3:tss+58][::-1]})
			tss_start_dict[tss] = array_of_features[tss].copy()

	return tss_start_dict



def write_csv(tss_dict_fwd,tss_dict_rev,output):
	tss_dict_fwd.update(tss_dict_rev)
	tss_dict = defaultdict(lambda:0)
	for tss in tss_dict_fwd:
		for key in tss_dict_fwd[tss]:
			tss_dict_fwd[tss][key] = str(tss_dict_fwd[tss][key]).translate(None,"['']")
		tss_dict[tss]= tss_dict_fwd[tss]
	df = pandas.DataFrame.from_dict(tss_dict,orient="index")
	pandas.DataFrame.to_csv(df,output)
	
def main(bamfile,genbankfile,output,window_size,step_height,min_p_value,min_slope):
	array_of_features_plus,array_of_features_minus,genome = parse_genbank(genbankfile) #Parsen des Genbank-File
	coverage_plus, coverage_minus = read_and_parse_bam_to_array(bamfile,genome) #Parsen des BAM-File und berechnen der Coverage-Arrays
	tss_array_plus = tss_detection_plus(coverage_plus,window_size,step_height,min_p_value,min_slope) #TSS-Detection
	tss_array_minus = tss_detection_minus(coverage_minus,window_size,step_height,min_p_value,min_slope) #TSS-Detection
	tss_dict_fwd = get_start_plus(array_of_features_plus,tss_array_plus,genome) #Unter Verwendung der Genomannotation, werden die detektierten TSS,wenn möglich, den zugehörigen Features zugewiesen.
	tss_dict_rev = get_start_minus(array_of_features_minus,tss_array_minus,genome)
	write_csv(tss_dict_fwd,tss_dict_rev,output)#Schreiben des Output-CSV-File

parser = argparse.ArgumentParser(description="Transcription-Start-Site Detection")
parser.add_argument('-b',action="store",dest="bamfile",type=str,required=True,help="Specify which BAM-File should be analyzed")
parser.add_argument('-g',action="store",dest="genbankfile",type=str,required=True,help="Specify .gb-File containing annotated features and reference genome sequence")
parser.add_argument('-o',action="store",dest="output",type=str,required=True,help="Specify Output-File (.csv) which contains detected TSS")
parser.add_argument('-w',action="store",dest="window_size",default=10,type=int,help="Specify regression window size")
parser.add_argument('-p',action="store",dest="min_p_value",default=0.0065,type=float,help="Specify regression p-value")
parser.add_argument('-s',action="store",dest="min_slope",default=0,type=int,help="Specify minimal regression slope")
parser.add_argument('-t',action="store",dest="step_height",default=20,type=int,help="Specify step height threshold")
args = parser.parse_args()

main(args.bamfile,args.genbankfile,args.output,args.window_size,args.step_height,args.min_p_value,args.min_slope)

	


	

