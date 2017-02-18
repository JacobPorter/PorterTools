#!/usr/bin/python
import optparse
import SeqIterator
import Constants
import sys

"""
A program for extracting sequences from a FASTA file
where the record ID matches a specified string.
This is for extracting primary assembly chromosomes
from a reference genome.
"""

def extractRef(in_file, match_string, out_file_d):
    ref_iterator = SeqIterator.SeqIterator(in_file, file_type = Constants.FASTA)
    ref_writer = SeqIterator.SeqWriter(out_file_d, file_type=Constants.FASTA)
    for fasta_record in ref_iterator:
        if match_string in fasta_record[0]:
            ref_writer.write(fasta_record)

def main():
    p = optparse.OptionParser()
    p.add_option('--fasta', '-f', help = 'The fasta file input.', default = None)
    p.add_option('--match', '-m', help = 'String to match on the chromosome identity', default = None)
    options, _ = p.parse_args()
    extractRef(options.fasta, options.match, sys.stdout)

if __name__ == '__main__':
    main()