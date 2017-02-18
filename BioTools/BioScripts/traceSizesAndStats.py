#!/usr/bin/python
from Bio import SeqIO

import math, optparse, sys

def main():
    
    p = optparse.OptionParser()
    
    p.add_option('--file', '-f' , help='The fasta or fastq file to check.')
    p.add_option('--bucket_size', '-b', default='10', help='The size of the buckets')
    p.add_option('--num_buckets', '-n', default='100', help='The number of buckets to have.  num_buckets * bucket_size = largest possible sequence size')
    p.add_option('--type', '-t', default='fasta', help='The type of file to count. Must be either "fasta" or "fastq"')
    p.add_option('-v', action='store_true', dest='validate', help='When turned on, checks if the sequences consist only of valid chars: ACTGNactgn')
    options, _ = p.parse_args()
    
    findHist(options.file, int(options.bucket_size), int(options.num_buckets), options.type, options.validate)
    

def findHist(filein, b_size, b_num, t, validate):
    seq_lens = []
    
    hist = [0] * b_num
    
    divisor = b_size
    
    valid_seq_char = "ACTGNactgn"
    valid_recs = True
    
    print "File %s being checked" % filein
    max_length = b_size * b_num
    print "Type: %s\tBucket size: %s\tNumber of buckets: %s\tMax length of sequence: %s" % (t, str(b_size), str(b_num), str(max_length))
    if validate:
        print "Validating sequence correctness"
    else:
        print "Validation turned off"
    print "---------------"
    
    sys.stdout.flush()
    count = 0
    for rec in SeqIO.parse(filein, t):
        if count % 1000000 == 0:
            print "%s sequences processed." % str(count)
            sys.stdout.flush()
        if validate:
            for l in rec.seq:
                if l not in valid_seq_char:
                    valid_recs = False
                    print "Invalid character found!!  " + str(rec.id)
                    sys.stdout.flush()
        len_seq = len(rec.seq)
        hist[int(len_seq / divisor)] += 1
        seq_lens.append(len_seq)
        count += 1
        
    seq_sum = sum(seq_lens)
    seq_len = len(seq_lens)
    seq_avg = seq_sum / (seq_len + 0.0)
    
    std_lst = map(lambda x: math.pow(x-seq_avg, 2), seq_lens)
    
    std_dev_seq = math.sqrt((sum(std_lst) / (seq_len + 0.0)))
    
    if not valid_recs and validate:
        print "Validation check FAILED.  Invalid chars found somewhere.  Please check."
    elif valid_recs and validate:
        print "Validation check PASSED."
        
    
    print "Total sequences found: " + str(count)
    #print "Smallest sequence size: " +  
    #print "Largest sequence size: " + 
    print "len " + str(seq_len)
    print "avg " + str(seq_avg)
    print "std dev " + str(std_dev_seq)
    print hist

    sys.stdout.flush()
if __name__ == '__main__':
    main()
