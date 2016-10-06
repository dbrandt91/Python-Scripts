# -*- coding: utf-8 -*-

from sys import argv
import pysam
import re
import argparse
import numpy as np

#BAM-File einlesen


def main(bamfile,mismatches):
    index_of_dot = bamfile.index(".bam")
    filename_good = bamfile[:index_of_dot]+"_good_mappings.bam"
    filename_bad = bamfile[:index_of_dot]+"_bad_mappings.bam"
    read_erased_count = 0
    read_kept_count = 0
    header = pysam.AlignmentFile(bamfile,"rb").header

    outfile_good = pysam.AlignmentFile(filename_good, "wb", header=header)
    outfile_bad = pysam.AlignmentFile(filename_bad, "wb", header=header)

    for read in pysam.AlignmentFile(bamfile,"rb"):
		read_stats = 0
		read_cigar = [i[1] for i in read.cigartuples if i[0] != 0]
		
		if read.get_tag("XM") > mismatches or len(read_cigar) > 0:
			if len(read_cigar) > 0: #indel criteria
				read_stats += read.get_tag("XM")
				for i in range(0,len(read_cigar)):
					if read_cigar[i] > mismatches or read_stats >= mismatches: #get number of indels 
						#print read.cigarstring, read.cigartuples,read.get_tag("XM")
						outfile_bad.write(read)
						read_erased_count += 1
						break
					elif read_cigar[i] < mismatches and read_stats < mismatches:
						read_stats += read_cigar[i]
				else:
					#print read.cigarstring, read.cigartuples,read.get_tag("XM")
					outfile_good.write(read)
					read_kept_count += 1	
						
			elif read.get_tag("XM") > mismatches: #mismatch criteria
				#print read.cigarstring, read.cigartuples,read.get_tag("XM")
				outfile_bad.write(read)
				read_erased_count += 1
		else:
			#print read.cigarstring, read.cigartuples,read.get_tag("XM")
			outfile_good.write(read)
			read_kept_count += 1
            
    print "Es wurden %d von %d Reads entfernt. Das entspricht %d %%" %(read_erased_count,read_erased_count+read_kept_count,(float(read_erased_count)/float(read_kept_count))*100)
    

parser = argparse.ArgumentParser(description="BAM-Quality-Filter for mapped reads")
parser.add_argument('-n',action="store",dest="mismatches",default=1,type=int,help="Specify the number of accepted mismatches")
parser.add_argument('-b',action="store",dest="bamfile",type=str,required=True,help="Specify which BAM-File should be analyzed")
parser.add_argument('-f',action="store",dest="firstbase",type=bool,required=False,help="Set Parameter = True if a mismatch in first base should lead to disposal of the respective read")
args = parser.parse_args()

main(args.bamfile,args.mismatches)
