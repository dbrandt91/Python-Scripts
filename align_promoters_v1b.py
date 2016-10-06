# -*- coding: utf-8 -*-

#Als Input dienen auschließlich die Improbizer-Results in einer Textdatei. Diese werden einfach komplett
#aus dem Browser kopiert.

from Bio import motifs,SeqIO
import pandas
import re
from sys import argv
script,txt= argv

list_of_seqs=[]
improbizer_result = []
index =[]
counter = 0
with open(txt,"r") as infile:
	for line in infile:
		if re.search("[A,C,G,T]",str(line)):
			line = str(line).split()[2]
			list_of_seqs.append(line.upper())
			line = re.sub("[a,c,g,t]","",line) #Improbizer-Results auf Motiv beschränken
			improbizer_result.append(line)
			index.append(list_of_seqs[counter].index(line)) #Motiv in Promotorseqs suchen
			counter += 1
		else:
			pass

#Ausrichten der Sequenzen
minimum = min(index)
#Basen am Anfang entfernen 
i = 0
for line in list_of_seqs:
	diff = index[i]-minimum
	if diff > 0:
		for k in range(diff):
			line = line[1:]
		list_of_seqs[i] = line
	else:
		pass
	i += 1
##Minimale String-Länge bestimmen
c= 1000
for line in list_of_seqs:
	if len(line) < c:
		c = len(line)
	else:
		pass
#Basen am Ende entfernen
i=0
for line in list_of_seqs:
	diff = len(line)-c
	if diff > 0:
		for k in range(diff):
			line= line[:-1]
		list_of_seqs[i] = line
	else:
		pass
	i += 1

with open("aligned_promoters.fasta","w") as outfile:
	for line in enumerate(list_of_seqs):
		outfile.write(">%d\n%s\n" %(line[0],line[1]))

#Motiv-Objekt berechnen
try:
	motif_object = motifs.create(list_of_seqs,alphabet=None)
except KeyError:
	next
#Motiv-Objekt in String umwandeln
motif_string = motif_object.format("pfm")
#Erstellen eines Weblogos
motif_object.weblogo("weblogo_test.eps",format="EPS")
#print motif_string
print motif_string
print motif_object.consensus
