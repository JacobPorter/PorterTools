#!/usr/bin/python
import optparse, math, sys

def main():
    p = optparse.OptionParser()
    
    p.add_option('--input', '-i' , help='The recovered input sequence file')
    p.add_option('--left', '-l', help='The left file for output')
    p.add_option('--right', '-r', help='The right file for output')
    p.add_option('--original', '-o', help='The original (recovered) sequence file.')
    p.add_option('--num', '-n', help='The number of mismatches to select')
    options, _ = p.parse_args()
    
    select(options.input, options.left, options.right, options.original, options.num)
    
    
def select(inputs, l, r, o, num):
    f_i = open(inputs, 'r')
    f_l = open(l, 'w')
    f_r = open(r, 'w')
    f_o = open(o, 'w')
    num = int(num)
    print "Searching for sequence pairs with %s mismatches.  Processing the input file %s and writing to the output files left: %s and right: %s" % (str(num), inputs, l, r)
    sys.stdout.flush()
    columns = ["id", "left", "aligned", "right", "percent", "num_missed", "misses"]
    count = 0
    for line in f_i:
        line_dict = dict(zip(columns, line.split("\t")))
        if int(line_dict["num_missed"]) == num:
            count += 1
            (seq1, seq2) = cleanse(line_dict["left"], line_dict["right"])
            ids = ">" + line_dict["id"] + "#"
            f_l.write(ids + "L" + "\n" + seq1 + "\n")
            f_r.write(ids + "R" + "\n" + seq2 + "\n")
            f_o.write(ids + "O" + "\n" + line_dict["aligned"].replace("-", "") + "\n")
            if count % 1000 == 0:
                f_l.flush()
                f_r.flush()
                f_o.flush()
    print "Number of sequences found: %s" % str(count)
    f_i.close()
    f_l.close()
    f_o.close()
    f_r.close()
            
        
def cleanse(seq1, seq2):
    len1 = len(seq1)
    len2 = len(seq2)
    if len1 != len2:
        diff = int(math.fabs(len1-len2))
        if len1 < len2:
            seq2 = seq2[0:len2-diff]
        else:
            seq1 = seq1[0:len1-diff]
    return (seq1, seq2)

if __name__ == '__main__':
    main()