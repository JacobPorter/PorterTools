#!/usr/bin/python
import SeqIterator
import optparse
import sys
import datetime

"""
This gets the first n records in a sequence record file and writes them to the file specified.
@author: Jacob Porter
"""

def getFirstRecords(file_in, output, record_type, number):
    input_iterator = SeqIterator.SeqIterator(file_in, file_type = record_type)
    output_writer = SeqIterator.SeqWriter(open(output, 'w'), file_type = record_type)
    counter = 0
    for record in input_iterator:
        if counter >= number:
            break
        output_writer.write(record)
        counter += 1
    output_writer.flush()


def main():
    now = datetime.datetime.now()
    p = optparse.OptionParser()
    p.add_option('--input','-i',help='Sequence file to partition')
    p.add_option('--output', '-o', help='File prefix for writing to.  File will be overwritten if it exists.')
    p.add_option('--num','-n',help='Number of sequences to get.')
    p.add_option('--type', '-t', help='The file type to partition.  Must be fasta, fastq, or sam', default='fastq')
    options,_ = p.parse_args()
    getFirstRecords(options.input, options.output, options.type, int(options.num))
    later = datetime.datetime.now()
    sys.stderr.write("Done!  The process started at %s and took %s time.\n" % (str(now), str(later - now)))
    
    
if __name__ == "__main__":
    main()