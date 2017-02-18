#!/usr/bin/python
import optparse
import datetime
import sys
import SeqIterator
import Constants
import IndentedHelpFormatterWithNL

"""
Script for removing fastq or fasta records with Ns in them.
"""

def removeDuplicates(reads, file_type = Constants.FASTQ):
    inputReads = SeqIterator.SeqIterator(reads, file_type = file_type)
    outputWriter = SeqIterator.SeqWriter(sys.stdout, file_type = file_type)
    counter_N = 0
    counter_all = 0
    duplicateDict = {}
    for record in inputReads:
        counter_all += 1
        if record[0] in duplicateDict:
            counter_N += 1
        else:
            duplicateDict[record[0]] = 0
            outputWriter.write(record)
            if counter_all <= 10:
                sys.stderr.write(record[1] + "\n")
                sys.stderr.flush()
                outputWriter.flush()
    outputWriter.flush()
    return (counter_all, counter_N)

def main():
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <reads>\n\nExample:\n %prog -t fasta ./reads.fa 1> reads.noN.fa 2> reads.err\n"
    version = "%prog " + "0.1.0"
    description = """
    This program removes all records that contain wild card 'N' characters.  
    It can process either fasta or fastq reads.  
    The records are printed to standard out, but program information is printed to standard error."""
    epilog = Constants.creation_string
    p = optparse.OptionParser(usage = usage, version = version, description = description, epilog = epilog, 
                              formatter = IndentedHelpFormatterWithNL.IndentedHelpFormatterWithNL())
    p.add_option('--type' , '-t', help= "The type of file to be processed.  This must be either 'fasta' or 'fastq' [default: %default]", default=Constants.FASTQ)
    options, args = p.parse_args()
    file_type = options.type.lower()
    if file_type != Constants.FASTA and file_type != Constants.FASTQ:
        p.error("The file type was incorrect.")
    sys.stderr.write("Removing duplicates from " + args[0] + "\n")
    sys.stderr.flush()
    (counter_all, counter_N) = removeDuplicates(args[0], file_type = file_type)
    later = datetime.datetime.now()
    sys.stderr.write("Total records: " + str(counter_all) + " Records removed: " + str(counter_N) + "\n")
    sys.stderr.write("Time to process the file: " + str(later - now) + "\n")
    sys.stderr.flush()

if __name__ == "__main__":
    main()