from Bio import SeqIO
#from Bio.SeqIO.QualityIO import PairedFastaQualIterator
#import os
#import string
#import itertools
#from os import listdir
#from os.path import isfile, join
import math

#start_path = "/home/jsporter/"

#outfile = open("SC_CHR_22_7340_CLEAN.fastq","w")


valid_seq_char = "ACTGNactgn-."
valid_recs = True

for rec in SeqIO.parse("/home/jsporter/Real_Data/chr22_hg19.fa", "fasta"):
    for l in rec.seq:
        if l not in valid_seq_char:
            valid_recs = False
            print "Invalid character found!!  " + str(rec.id)
    

if not valid_recs:
    print "Invalid chars found somewhere.  Please check."

