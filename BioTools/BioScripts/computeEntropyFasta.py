#!/usr/bin/python
from Bio import SeqIO

import math, optparse, sys

def main():
    
    p = optparse.OptionParser()
    
    p.add_option('--file', '-f' , help='The fasta or fastq file to check.')
    p.add_option('--type', '-t', default='fasta', help='The type of file to count. Must be either "fasta" or "fastq"')
    p.add_option('-i', action='store_true', dest='ignore_ambig', help="When turned on, ignores sequences with ambiguous positions.")
    p.add_option('-v', action='store_true', dest='verbose', help="Prints out all frequency and entropy information for each sequence.")
    
    options, _ = p.parse_args()
    
    findEntropy(options.file, options.type, options.ignore_ambig, options.verbose)
    

def findEntropy(filein, t, ignore_ambig, verbose):
    bucket_num = 20
    bucket_size = 10
    print "Doing entropy calculations for the " + t + " file: " + filein
    if verbose:
        print "All frequency information is in the order: ACTGN"
    
    seq_count = 0
    total_seq = 0
    global_freq = [0, 0, 0, 0, 0]
    global_len = 0
    entropy_totals = []
    lookup = {"A":0, "C":1, "T":2, "G":3, "N":4}
    print ""
    sys.stdout.flush()
    my_hist = [0] * (bucket_num + 1)
    for rec in SeqIO.parse(filein, t):
        total_seq += 1
        local_freq = [0, 0, 0, 0, 0]
        local_len = 0
        for l in rec.seq:
            local_freq[lookup[l.upper()]] += 1
            local_len += 1
        if ignore_ambig and local_freq[lookup["N"]] > 0:
            continue
        else:
            local_entropy = 0
            for key in lookup:
                global_freq[lookup[key]] += local_freq[lookup[key]]
                if local_freq[lookup[key]] != 0:
                    my_freq = local_freq[lookup[key]] / (0.0 + local_len)
                    local_entropy += my_freq * math.log(my_freq, 2)
            local_entropy = local_entropy * -1
            try:
                my_hist[int(local_entropy * bucket_size)] += 1
            except IndexError:
                print "Record %s not added to histogram.  Entropy out of range" % str(rec.id)
            entropy_totals.append(local_entropy)
            global_len += local_len 
            seq_count += 1
            if verbose:
                print "Record ID: %s\tFreq: %s\tLength: %s\tSequence Entropy: %s" % (str(rec.id), str(local_freq), str(local_len), str(local_entropy))
                sys.stdout.flush()
    
    print ""
    entropy_sum = sum(entropy_totals)
    if verbose:
        print "Entropy sum: " + str(entropy_sum)
    entropy_avg = entropy_sum / (seq_count + 0.0)
    
    std_lst = map(lambda x: math.pow(x-entropy_avg, 2), entropy_totals)
    
    std_dev_entropy = math.sqrt((sum(std_lst) / (seq_count + 0.0)))
    print "Total number of nucleotides for sequences where entropy was calculated: " + str(global_len)
    print "Global count of nucleotides: " + str(global_freq)
    print "Total sequences used in calculation: " + str(seq_count)
    print "Total sequences ignored: " + str(total_seq - seq_count)
    sys.stdout.flush()
    global_entropy = 0.0
    for key in lookup:
        if local_freq[lookup[key]] != 0:
            my_freq = global_freq[lookup[key]] / (0.0 + global_len)
            global_entropy += my_freq * math.log(my_freq, 2)
    global_entropy = global_entropy * -1
    print "Entropy for all used sequences: " + str(global_entropy)
    print "Average sequence entropy: " + str(entropy_avg)
    sys.stdout.flush()
    print "Entropy Max: " + str(max(entropy_totals))
    sys.stdout.flush()
    print "Entropy Min: " + str(min(entropy_totals))
    sys.stdout.flush()
    print "Entropy Standard Deviation: " + str(std_dev_entropy)
    sys.stdout.flush()
    print ""
    print "Entropy Histogram: "
    print my_hist
    sys.stdout.flush()
    
    
    
if __name__ == '__main__':
    main()
