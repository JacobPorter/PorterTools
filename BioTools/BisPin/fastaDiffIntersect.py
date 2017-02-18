#!/usr/bin/python
import sys, optparse, re

def main():
    p = optparse.OptionParser()
    p.add_option('--file1', '-1' , help='The first fasta file to compare')
    p.add_option('--file2', '-2', help='The second fasta file to compare')
    p.add_option('--diff', '-d', default='diff.fa', help='The file prefix to store differences to.  Differences for file 1 and file 2 will be stored separately denoted by a 1 and 2 suffix')
    p.add_option('--intersect', '-i', default= 'inter.fa', help='The file to store intersections to')
    p.add_option('--reg_ex', '-r', help='The regular expression used to match sequence identifiers.  Identifiers are stripped of leading and trailing white space before comparison.  This feature be used to get rid of extra characters.', default='(.*)')
    p.add_option('-s', action='store_true', dest='sequences', help="When turned on, check the identity of fasta entries based on the sequences rather than the identifier that starts with '>'.")
    p.add_option('-c', action='store_true', dest='count_only', help="When turned on, only counts the intersections and differences and does not save them to a fasta file.")
    options, _ = p.parse_args()
    findDiffsIntersects(options.file1, options.file2, options.sequences, options.diff, options.intersect, options.count_only, options.reg_ex)
    
def writeFile(fd, identity, seq):
    fd.write(identity + "\n")
    fd.write(seq + "\n")
    fd.flush()
    
def findDiffsIntersects(f1, f2, seq, diff, inter, count_only, reg_ex):
    d1 = {}
    #Use the regular expression '(.+)#0/' for hairpin data
    my_compare = re.compile(reg_ex)
    #Process first file
    if seq:
        seq_info = "Using sequence information for fasta entry identity."
    else:
        seq_info = "Using fasta identifier for fasta entry identity."
    sys.stdout.write(seq_info + "\nUsing first fasta file for dictionary creation: " + f1 + "\n")
    sys.stdout.flush()
    f1_fd = open(f1)
    identity = ""
    my_seq = ""
    f1_count = 0
    for line in f1_fd:
        if line[0] == ">" and identity != "":
            if seq:
                d1[my_seq] = identity
            else:
                d1[identity] = my_seq
            identity = my_compare.match(line.strip()).groups()[0]
            my_seq = ""
            f1_count += 1
        elif line[0] == ">":
            identity = my_compare.match(line.strip()).groups()[0]
        else:
            my_seq += line.strip()
    if identity != "":
        if seq:
            d1[my_seq] = identity
        else:
            d1[identity] = my_seq
        f1_count += 1
    f1_fd.close()
    sys.stdout.write("Sequences found in the first file: " + str(f1_count) + "\nProcessing second file for differences and intersections: " + f2 + "\n")
    sys.stdout.flush()
    #Process second file
    f2_fd = open(f2)
    if not count_only:
        diff1_fd = open(diff + ".1.diff", "w")
        diff2_fd = open(diff + ".2.diff", "w")
        inter_fd = open(inter, "w")
    diffs_count_1 = 0
    diffs_count_2 = 0
    inter_count = 0
    identity = ""
    my_seq = ""
    f2_count = 0
    for line in f2_fd:
        if line[0] == ">" and identity != "":
            if seq:
                if my_seq in d1:
                    if not count_only:
                        writeFile(inter_fd, d1[my_seq] + " " + identity, my_seq)
                    inter_count += 1
                    d1.pop(my_seq, None)
                else:
                    if not count_only:
                        writeFile(diff2_fd, identity, my_seq)
                    diffs_count_2 += 1
            else:
                if identity in d1:
                    if not count_only:
                        writeFile(inter_fd, identity, my_seq)
                    inter_count += 1
                    d1.pop(identity, None)
                else:
                    if not count_only:
                        writeFile(diff2_fd, identity, my_seq)
                    diffs_count_2 += 1
            identity = my_compare.match(line.strip()).groups()[0]
            my_seq = ""
            f2_count += 1
        elif line[0] == ">":
            identity = my_compare.match(line.strip()).groups()[0]
        else:
            my_seq += line.strip()
    if identity != "":
        if seq:
            if my_seq in d1:
                if not count_only:
                    writeFile(inter_fd, d1[my_seq] + " " + identity, my_seq)
                inter_count += 1
                d1.pop(my_seq, None)
            else:
                if not count_only:
                    writeFile(diff2_fd, identity, my_seq)
                diffs_count_2 += 1
        else:
            if identity in d1:
                if not count_only:
                    writeFile(inter_fd, identity, my_seq)
                inter_count += 1
                d1.pop(identity, None)
            else:
                if not count_only:
                    writeFile(diff2_fd, identity, my_seq)
                diffs_count_2 += 1
        f2_count += 1
    f2_fd.close()
    if not count_only:
        inter_fd.close()
    sys.stdout.write("Processed " + str(f2_count) + " sequences for the second file.\n")
    sys.stdout.flush()
    #Process remaining sequences from the first file
    for key in d1:
        if not count_only:
            if seq:
                writeFile(diff1_fd, d1[key] , key)
            else:
                writeFile(diff1_fd, key, d1[key])
        diffs_count_1 += 1
    if not count_only:
        diff1_fd.close()
        diff2_fd.close()
    sys.stdout.write("Found " + str(inter_count) + " intersecting fasta entries \nFound " + str(diffs_count_1) + " different fasta entries in first file " + f1 + "\nFound " + str(diffs_count_2) + " different fasta entries in second file " + f2 + "\n")
    sys.stdout.flush()

if __name__ == '__main__':
    main()