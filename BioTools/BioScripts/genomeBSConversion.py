#!/usr/bin/python
from Bio import SeqIO
from Bio.Seq import MutableSeq
import optparse, sys, os


def main():
    
    p = optparse.OptionParser()
    
    p.add_option('--in_file', '-i' , help='The file to convert.')
    p.add_option('--out_file', '-o', help='The prefix of the file to save.')
    p.add_option('--type', '-t', default='fasta', help='The type of file: fasta or fastq.  Default is fasta')
    options, _ = p.parse_args()
    
    convert_wrapper(options.in_file, options.out_file, options.type)

def convert_wrapper(filein, fileout, t):
    lettersCT = (('C','c'), ('T','t'))
    convert(filein, fileout, t, lettersCT)
    lettersGA = (('G','g'), ('A','a'))
    convert(filein, fileout, t, lettersGA)
    
   
def convert(filein, fileout, t, letters):
    suffix_str = letters[0][0] + letters[1][0] + "." + t
    print "Computing %s conversion for file: %s" % (suffix_str, filein)
    sys.stdout.flush()
    my_recs = []
    
    count = 0
    for rec in SeqIO.parse(filein, t):
        print "Processing record %s." % str(count)
        sys.stdout.flush()
        my_seq = rec.seq.tomutable()
        for i in xrange(len(my_seq)):
            l = my_seq[i]
            if l in letters[0]:
                my_seq[i] = letters[1][letters[0].index(l)]
        rec.seq = my_seq.toseq()
        my_recs.append(rec)
    
    output_handle = open(os.path.join(fileout + suffix_str), "w")
    SeqIO.write(my_recs, output_handle, t)
    output_handle.close()
    print "Saving file to %s" % (fileout + suffix_str)
    sys.stdout.flush()
    count += 1
    
if __name__ == '__main__':
    main()
