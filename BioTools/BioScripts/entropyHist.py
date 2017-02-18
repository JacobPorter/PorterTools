#!/usr/bin/python
from Bio import SeqIO
import cPickle as Pickle
import optparse, sys


def main():
    p = optparse.OptionParser()
    p.add_option('--fasta', '-f', help='Fasta file that has sequence ids')
    p.add_option('--dict', '-d', help='Dictionary file from createEntDict that has sequence ids and entropy info')
    p.add_option('--bucket_num', '', default='20', help='The number of buckets in the histogram.')
    p.add_option('--bucket_size', '', default='10', help='The number of buckets between 0 and 1.')
    
    options, _ = p.parse_args()
    
    findHist(options.fasta, options.dict, int(options.bucket_num), float(options.bucket_size))
    
def findHist(fasta_file, dict_file, bucket_num, bucket_size):
    print "#Using fasta file: " + fasta_file
    print "#Using dictionary file: " + dict_file
    sys.stdout.flush()
    #fasta_fd = open(fasta_file)
    ent_dict = Pickle.load(open(dict_file))
#     for key in ent_dict:
#         print key
#         sys.stdout.flush()
    my_hist = [0] * (bucket_num + 1)
    count = 0
    print ""
    for rec in SeqIO.parse(fasta_file, "fasta"):#for line in fasta_fd:
        if "N" in str(rec.seq).upper():
            print "#Ambiguous record found: " + str(rec.id)
            sys.stdout.flush()
            continue
        else:
            count += 1
        #rec.id
        #if line[0] == '>':
            
            #f_id = line[1:].strip()
            entropy = int(float(ent_dict[str(rec.id)][2].split(":")[1].strip()) * bucket_size)
            try:
                my_hist[entropy] += 1
            except IndexError:
                print "%s produced index error with entropy: %s" % (str(rec.id), str(entropy))
                sys.stdout.flush()
                
    
    print "#Records processed: " + str(count)
    print ""
    sys.stdout.flush()
    print "#Entropy histogram:"
    print ""
    sys.stdout.flush()
    print my_hist
    sys.stdout.flush()
    #fasta_fd.close()
    
if __name__ == '__main__':
    main()
    