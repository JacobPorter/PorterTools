#!/usr/bin/python
from Bio import SeqIO

import math, optparse, sys

def main():
    
    p = optparse.OptionParser()
    
    p.add_option('--file_in', '-i' , help='The fasta or fastq file to extract from.')
    p.add_option('--file_out', '-o' , help='The file to save the extracted sequences to.')
    p.add_option('--type', '-t', default='fasta', help='The type of file.  Must be either "fasta" or "fastq"')
    p.add_option('--length', '-l', help='The sequence length of extracted sequences.')
    options, _ = p.parse_args()
    
    extractSeq(options.file_in, options.file_out, options.type, int(options.length))
    

def extractSeq(filein, fileout, t, length):
    print "Extracting sequences from %s file from %s " % (t, filein)
    print "Sequence size to be extracted: %s" % str(length)
    print "Saving sequences to %s" % fileout
    print "---------------"
    
    sys.stdout.flush()
    input_handle = open(filein, "rU")
    output_handle = open(fileout, "w")
    
    count = 0
    write_count = 0
    for rec in SeqIO.parse(input_handle, t):
        if count % 1000000 == 0:
            print "%s sequences processed." % str(count)
            sys.stdout.flush()
        len_seq = len(rec.seq)
        if len_seq == length:
            write_count += 1
            SeqIO.write(rec, output_handle, t)
        count += 1
    
    print "Number of sequences in input file: %s" % count
    print "Number of sequences extracted: %s" % write_count
    sys.stdout.flush()
    input_handle.close()
    output_handle.close()
        
if __name__ == '__main__':
    main()
