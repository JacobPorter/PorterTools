#!/usr/bin/python
import SeqGenerator, optparse, sys, random



def main():
    p = optparse.OptionParser()
    p.add_option("--in_file", "-i", help="Input file.  Either fasta or fastq file.")
    p.add_option("--out", "-o", help="Output file to write to.  Will be appended if it exists.")
    p.add_option("--type", "-t", default="fastq", help="Type of input file.  Must be either 'fasta' or 'fastq.'")
    p.add_option("--length", "-l", default="100", help="Length of each sequence to sample.")
    p.add_option("--amount", "-a", default="100000", help="")
    options, _ = p.parse_args()
    start(options.in_file, options.out, options.type, int(options.length), int(options.amount))
    print "Done!"

def start(in_f, out_f, t, length, amount):
    print "Sampling from %s and writing to %s.\nLength %s.  Amount %s.  Type %s" % (in_f, out_f, str(length), str(amount), t)
    sys.stdout.flush()
    (seq_lst, seq_count, nt_count) = getLengths(in_f, t, length)
    print seq_lst
    print "Sequence Count", seq_count, "Possible positions to sample from", nt_count
    sys.stdout.flush()
    rands = [random.randint(0, nt_count-1) for i in xrange(amount) ]
    rands.sort()
#     if debug():
#         print rands
    getRandomSequences(rands, in_f, out_f, t, length, amount, seq_lst, seq_count, nt_count)

def getLengths(file_str, t, length):
    seq_lst = []
    seq_count = 0
    nt_count = 0
    seqs_gen = SeqGenerator.SeqGenerator(file_str, file_type=t)
    for rec in seqs_gen:
        seq_len = len(rec[1])
        if seq_len < length:
            continue
#         if debug():
#             rec_tpl = (seq_count, rec[0], seq_len, nt_count, nt_count + seq_len - length + 1)
#             seq_lst.append(rec_tpl)
        nt_count += (seq_len - length + 1)
        seq_count += 1
    return (seq_lst, seq_count, nt_count)

def getID(start_id, num, b, e):
    return start_id + "_" + str(num) + "_b" + str(b) + "_e" + str(e)
     
def getRandomSequences(rands, in_f, out_f, t, length, amount, seq_lst, seq_count, nt_count):
    my_nt_count = 0
    my_rand_count = 0
    seqs_gen = SeqGenerator.SeqGenerator(in_f, file_type=t)
    fd_seq_write = open(out_f, 'w')
    seq_writer = SeqGenerator.SeqWriter(fd_seq_write, file_type=t)
    for rec in seqs_gen:
        seq_len = len(rec[1])
        if seq_len < length:
            continue
#         if debug():
#             print "my_nt_count", my_nt_count, "my_nt_count and length", my_nt_count + seq_len - length + 1, "next random position", rands[my_rand_count], "continue", str(my_nt_count + seq_len - length + 1 < rands[my_rand_count])
        if my_nt_count + seq_len - length + 1 < rands[my_rand_count]:
            my_nt_count += seq_len - length + 1
            continue
        num = 0
        while my_nt_count + seq_len - length >= rands[my_rand_count]:
            b = rands[my_rand_count] - my_nt_count
            e = b + length
#             if debug():
#                 print "id", str(rec[0]), "begin", str(b), "end", str(e), "my_rand_pos", rands[my_rand_count], "my_nt_count", my_nt_count
            if t == "fastq":
                nextSeq = (getID(rec[0], num, b, e) , rec[1][b:e], rec[2][b:e] )
            elif t == "fasta":
                nextSeq = (getID(rec[0], num, b, e), rec[1][b:e])
            seq_writer.write(nextSeq)
            my_rand_count += 1
            num += 1
            if my_rand_count >= len(rands):
                return
        my_nt_count += seq_len - length + 1
    
    
if __name__ == "__main__":
    main()