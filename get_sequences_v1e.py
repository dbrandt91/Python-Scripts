#get_sequences.py
#extracts multiple sequences from fasta file by coordinates (Locus, Start, Stop, Strand) and saves them into seq_output.fasta
#Usage: python get_sequences.py <genome.fasta> <input.csv>
#Columns in CSV have to be named as: "Locus", "Start", "Stop" and "Strand [Fwd or Rev]"

import pandas
from Bio import SeqIO
from sys import argv
script,fasta,csv = argv

def get_sequences(fasta,csv):
	index = str(csv).index(".csv")
	output_file = str(csv)[0:index]+".fasta"
	df = pandas.read_csv(csv,sep=",")
	#df.columns = [c.replace(" ","_") for c in df.columns]
	starts = df.Start.tolist()
	stops = df.Stop.tolist()
	strand = df.Strand.tolist()
	locus = df.Locus.tolist()
	locus = map(lambda x: str(x).replace(" ",""),locus) #Leerzeichen in der Locus-Bezeichnung entfernen
	seqs = [None]*len(starts)
	output_txt = open(output_file,"w")
	with open(fasta,"rU") as fasta:
		for record in SeqIO.parse(fasta,"fasta"):
			genome_length = len(str(record.seq))
			for i in range(0,len(starts)):
				if locus[i] == "nan":
					next
				elif strand[i] == "Fwd":
					if int(starts[i]) > int(stops[i]):
						seqs[i] = str(record.seq[int(starts[i])-1:genome_length]).upper()+str(record.seq[0:int(stops[i])]).upper()
						print int(starts[i])-1,genome_length,int(stops[i])
					else:
						seqs[i] = str(record.seq[int(starts[i])-1:int(stops[i])]).upper()
				elif strand[i] =="Rev":
					if starts[i] < stops[i]:
						seqs[i] = str(record.seq[int(starts[i])-1:int(stops[i])].reverse_complement()).upper()
					elif stops[i] < starts[i]:
						seqs[i] = str(record.seq[int(stops[i])-1:int(starts[i])].reverse_complement()).upper()
				if seqs[i] == None:
					next 
				elif i == 0:
					output_txt.write(">%s\n%s" %(locus[i],seqs[i]))
				else:
					output_txt.write("\n>%s\n%s" %(locus[i],seqs[i]))
	output_txt.close()

get_sequences(fasta,csv)
