#!/usr/bin/python
import SeqGenerator
import optparse, sys

"""
Sorts the fasta file by the identifier.
"""

def main():
    p = optparse.OptionParser()
    p.add_option('--input','-i',help='Fasta file to sort')
    p.add_option('--output', '-o', help='File to write to.  File will be overwritten if it exists.')
    options,args = p.parse_args()
#    if len(args) == 0:
#        print "Use --help for help documentation"
#        return
    fasta_sort(options.input, options.output)
    print "Done!"

def fasta_sort(file_in, file_out):
    seq_g = SeqGenerator.SeqGenerator(file_in)
    print "Sorting sequences by id from file %s" % file_in
    sys.stdout.flush()
    seq_lst = [s for s in seq_g]
    seq_lst.sort(key=lambda seq: seq[0])
    print "Writing %d sorted sequences to file %s" % (len(seq_lst), file_out)
    sys.stdout.flush()
    seq_w = SeqGenerator.SeqWriter(open(file_out, 'w'))
    count = 0
    for s in seq_lst:
        count += 1
        seq_w.write(s)
        if count % 10000 == 0:
            seq_w.flush()
    print "Sequences written: %d" % count
    

if __name__ == '__main__':
    main()