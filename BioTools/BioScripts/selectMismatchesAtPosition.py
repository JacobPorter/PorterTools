#!/usr/bin/python
import optparse, math, sys

def main():
    p = optparse.OptionParser()
    
    p.add_option('--input', '-i' , help='The recovered input sequence file')
    p.add_option('--left', '-l', help='The left file for output')
    p.add_option('--right', '-r', help='The right file for output')
    p.add_option('--original', '-o', help='The original (recovered) sequence file.')
    p.add_option('--num', '-n', help='The number of mismatches to select')
    p.add_option('--begin', '-b', help='The beginning window for mismatches')
    p.add_option('--end', '-e', help='The end of the window where there are mismatches')
    options, _ = p.parse_args()
    
    select(options.input, options.left, options.right, options.original, options.num, options.begin, options.end)
    
def find_mismatch_locations(mismatch_string):
    end = len(mismatch_string)
    mm = mismatch_string[0:end-2].split(":")
    mm_lst = []
    for i in mm:
        mm_lst.append(i.split("_"))
    mm_lst.sort()
    #print mm_lst
    return mm_lst
    
def inside(my_misses, begin, end):
    ret_bool = True
    for i in my_misses:
        if int(i[0]) < begin or int(i[0]) > end:
            ret_bool = False
    return ret_bool


def outside(my_misses, begin, end):
    ret_bool = True
    for i in my_misses:
        if int(i[0]) >= begin and int(i[0]) <= end:
            ret_bool = False
    return ret_bool
    
def select(inputs, l, r, o, num, begin, end):
    f_i = open(inputs, 'r')
    f_l = open(l, 'w')
    f_r = open(r, 'w')
    f_o = open(o, 'w')
    f_l_o = open(l+"out", 'w')
    f_r_o = open(r+"out", 'w')
    f_o_o = open(o+"out", 'w')
    num = int(num)
    begin = int(begin)
    end = int(end)
    print "Searching for sequence pairs with %s mismatches.  Processing the input file %s and writing to the output files left: %s and right: %s" % (str(num), inputs, l, r)
    sys.stdout.flush()
    columns = ["id", "left", "aligned", "right", "percent", "num_missed", "misses"]
    count = 0
    for line in f_i:
        line_dict = dict(zip(columns, line.split("\t")))
        if int(line_dict["num_missed"]) == num:
            count += 1
            print line_dict["id"], line_dict["aligned"]
            my_misses = find_mismatch_locations(line_dict["misses"])
            print my_misses
            if inside(my_misses, begin, end): #int(my_misses[0][0]) >= begin and int(my_misses[-1][0]) <= end:
                print "inside"
                (seq1, seq2) = cleanse(line_dict["left"], line_dict["right"])
                ids = ">" + line_dict["id"] + "#"
                f_l.write(ids + "L_W" + "\n" + seq1 + "\n")
                f_r.write(ids + "R_W" + "\n" + seq2 + "\n")
                f_o.write(ids + "O_W" + "\n" + line_dict["aligned"].replace("-", "") + "\n")
            elif outside(my_misses, begin, end):
                print "outside"
                (seq1, seq2) = cleanse(line_dict["left"], line_dict["right"])
                ids = ">" + line_dict["id"] + "#"
                f_l_o.write(ids + "L_O" + "\n" + seq1 + "\n")
                f_r_o.write(ids + "R_O" + "\n" + seq2 + "\n")
                f_o_o.write(ids + "O_O" + "\n" + line_dict["aligned"].replace("-", "") + "\n")
            else:
                print "neither"
            sys.stdout.flush()
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