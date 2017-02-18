#!/usr/bin/python
import optparse
import SeqIterator
import Constants
import sys
import datetime

"""
"This program converts a FASTA file into a multiline format with the line length specified as an argument."
@author: Jacob Porter
@since: 01/14/2017
"""

logstr = "multiLineFasta: "

def outputFasta(input_file, output_file, lineLength):
    input_iterator = SeqIterator.SeqIterator(input_file)
    output_iterator = SeqIterator.SeqWriter(open(output_file, 'w'), file_type = Constants.FASTA, 
                                            line_toggle = True, line_length = lineLength)
    for record in input_iterator:
        output_iterator.write(record)
    sys.stderr.write("%sNumber of records processed:\t%s\n" % (logstr, str(input_iterator.records_processed())))

def main():
    usage = "usage: %prog [options] <reference_genome_file> <output_file_name>"
    description = "This program converts a fasta file into a multiline format with the line length specified as an argument."
    p = optparse.OptionParser(usage = usage, description = description)
    p.add_option('--lineLength','-l',help= 'Determines the length of each line', default = "60")
    options, args = p.parse_args()
    if len(args) != 2:
        p.error("A fasta file input and an output file name must be specified as the only required arguments.")
    try:
        lineLength = int(options.lineLength)
        if lineLength <= 0:
            raise ValueError
    except ValueError:
        p.error("The line length argument must be a positive integer.")
    sys.stderr.write("%sStarting multiLineFasta on input %s and output %s with line length %s on %s\n" % 
                     (logstr, args[0], args[1], str(lineLength), str(datetime.datetime.now())))
    outputFasta(args[0], args[1], lineLength)

if __name__ == "__main__":
    main()