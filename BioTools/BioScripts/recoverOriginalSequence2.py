#!/usr/bin/python
import SeqGenerator
import optparse
import math, sys, itertools #,nw


def main():
    p = optparse.OptionParser()
    
    p.add_option('--CtoT', '-c' , help='The fasta file with the C to T converted sequences.')
    p.add_option('--GtoA', '-g', help='The fasta file with the corresponding G to A converted sequences.')
    p.add_option('--outFile', '-o', help='The location and prefix for the output file to write to')
    options, _ = p.parse_args()
    
    recover(options.CtoT, options.GtoA, options.outFile)
    print "Done!"
    

def recover(CtoT_file, GtoA_file, outFile):
    out_fd = open(outFile, 'w')
    CtoT_Seqs = SeqGenerator.SeqGenerator(CtoT_file)
    GtoA_Seqs = SeqGenerator.SeqGenerator(GtoA_file)
    print "Aligning sequences from C-to-T file %s and G-to-A file %s and writing output to %s" % (CtoT_file, GtoA_file, outFile)
    sys.stdout.flush()
    count = 0
    for i,j in itertools.izip(CtoT_Seqs, GtoA_Seqs):
        if sameSeq(i,j):
            if count % 1000 == 0:
                print "Count: " + str(count)
                sys.stdout.flush()
            #print i[1]
            #print j[1]
            (my_align, _, _, percent_match) = align(i, j, out_fd)
            #print my_align, percent_match
            count += 1
        else:
            print "Error!!! Sequences do not have the same identifier.  %s, %s" % (i[0], j[0])
            sys.stdout.flush()
    print "Total sequences processed: %d" % count
    
#    for i in CtoT_Seqs:
#        GtoA_Seqs.reset()
#        #print str(i[0].split('#')[0])
#        for j in GtoA_Seqs:
#            if sameSeq(i,j):
#                if count % 1000 == 0:
#                    print count
#                    sys.stdout.flush()
#                #print i[1]
#                #print j[1]
#                (my_align, _, _, percent_match) = align(i, j, out_fd)
#                #print my_align, percent_match
#                count += 1
#                break

def align(i, j, out_fd):
    len1 = len(i[1])
    len2 = len(j[1])
    my_align = ""
    mismatches = ""
    mismatch_cnt = 0
    for pos in range(min(len1, len2)):
        ib = i[1][pos]
        jb = j[1][pos]
        if ib == jb:
            my_align += ib
        elif ib == "T" and jb == "C":
            my_align += "C"
        elif ib == "G" and jb == "A":
            my_align += "G"
        elif ib != jb:
            my_align += "*"
            mismatches += str(pos) + "_" + ib + jb + ":"
            mismatch_cnt += 1
    for pos in range(int(math.fabs(len1-len2))):
        my_align += "-"
    percent_match = (max(len1, len2) - mismatch_cnt - math.fabs(len1-len2)) / (max(len1, len2) + 0.0)
    out_fd.write(str(i[0].split('#')[0]) + "\t" +  str(i[1]) + "\t"  + my_align + "\t" + str(j[1]) + "\t" + str(percent_match) + "\t" + str(mismatch_cnt) + "\t" + str(mismatches) + "\n")
    out_fd.flush()
    return (my_align, mismatches, mismatch_cnt, percent_match)
        
   
#def align(i, j, out_fd):
#    (penalty, align1, align2) = nw.needleman_wunsch(-3, -1, i[1], j[1])
#    print align1
#    print align2
#    out_fd.write(str(i[0].split('#')[0]) + "\t" + str(penalty) + "\t" + str(i[1]) + "\t" + str(align1) + "\t" + str(j[1]) + "\t" + str(align2) + "\n")
#    out_fd.flush()
    
def sameSeq(i,j):
    return i[0].split('#')[0] == j[0].split('#')[0]

if __name__ == '__main__':
    main()