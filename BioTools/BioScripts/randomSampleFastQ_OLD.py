#!/usr/bin/python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import optparse, sys, gzip, random


def main():
    p = optparse.OptionParser()
    p.add_option("--in_file", "-i", help="Input file.  File must be gzipped  Either gzipped fasta or fastq file.")
    p.add_option("--out", "-o", help="Output file to write to.  Will be appended if it exists.")
    p.add_option("--type", "-t", default="fastq", help="Type of input file.  Must be either 'fasta' or 'fastq.'")
    p.add_option("--length", "-l", default="100", help="Length of each sequence to sample.")
    p.add_option("--amount", "-a", default="100000", help="")
    
    options, _ = p.parse_args()
    start(options.in_file, options.out, options.type, int(options.length), int(options.amount))

def start(in_f, out_f, t, length, amount):
    print ""
    sys.stdout.flush()
    (seq_lst, seq_count, nt_count) = getLengths(in_f, t, length)
    rands = [random.randint(0, nt_count-1) for i in xrange(amount) ]
    rands.sort()
    getRandomSequences(rands, in_f, out_f, t, length, amount, seq_lst, seq_count, nt_count)


def getLengths(file_str, t, length):
    seq_lst = []
    #(position, id, length)
    seq_count = 0
    nt_count = 0
    handle = gzip.open(file_str)
    for rec in SeqIO.parse(handle, t):
        seq_len = len(rec.seq)
        if seq_len < length:
            continue
        rec_tpl = (seq_count, rec.id, seq_len, nt_count, nt_count + seq_len - length + 1)
        seq_lst.append(rec_tpl)
        nt_count += (seq_len - length + 1)
        seq_count += 1
    #handle.gzip.close()
        
    return (seq_lst, seq_count, nt_count)

# def generateRandomPositions(nt_count, amount):
#     pass
     
def getRandomSequences(rands, in_f, out_f, t, length, amount, seq_lst, seq_count, nt_count):
    in_handle = gzip.open(in_f)
    out_handle = gzip.open(out_f, "w+")
    my_seq_count = 0
    my_nt_count = 0
    my_rand_count = 0
    
    for rec in SeqIO.parse(in_handle, t):
        seq_len = len(rec.seq)
        if seq_len < length:
            continue
        temp = my_seq_count
        num = 0
        while my_seq_count == temp:
            b = rands[my_rand_count] - my_nt_count
            e = b + length
            nextSeq = SeqRecord(rec.seq[b:e], id=rec.id + "__" + str(num))
            SeqIO.write(nextSeq, out_handle, t)
            my_rand_count += 1
            num += 1
            if (rands[my_rand_count] >= my_nt_count + seq_len - length + 1):
                my_nt_count += seq_len - length + 1
                my_seq_count += 1
    
    

if __name__ == "__main__":
    main()
