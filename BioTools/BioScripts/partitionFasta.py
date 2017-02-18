#!/usr/bin/python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys, optparse

def main():
    p = optparse.OptionParser()
    
    p.add_option('--input', '-i' , help='Input fasta file')
    p.add_option('--num', '-n', help='Number of partitions')
    options, _ = p.parse_args()
    
    partition(options.input, int(options.num))

def partition(filein, num):
    
    in_handle = open(filein, 'r')
    count = 0
    for rec in SeqIO.parse(in_handle, 'fasta'):
        count += 1
    in_handle.close()
    in_handle = open(filein, 'r')
    part_size = count / (num)
    part_num = 0
    part_count = 0
    part_file = open(filein + str(part_num), 'w')
    print "Working on partition file %s" % (str(part_num))
    sys.stdout.flush()
    for rec in SeqIO.parse(in_handle, 'fasta'):
        if part_count < part_size or part_num == num - 1:
            part_count += 1
            SeqIO.write([rec], part_file, 'fasta')
            # write rec to file
        else:
            part_count = 0
            part_num += 1
            part_file.close()
            part_file = open(filein + str(part_num), 'w')
            print "Working on partition file %s" % (str(part_num))
            sys.stdout.flush()
            SeqIO.write([rec], part_file, 'fasta')
            part_count += 1
    part_file.close()
    print "Done"
        
if __name__ == '__main__':
    main()