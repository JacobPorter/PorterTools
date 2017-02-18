#!/usr/bin/env python
 
# Needleman-Wunsch algorithm for calculating 
# optimal alignment of two DNA sequences. 
# http://evandempsey.wordpress.com/2013/01/08/needleman-wunsch-algorithm-for-dna-sequence-alignment/
# evandempsey
 
import sys
 
def read():
    """
    Read sequences from file.
    """
    f = open(sys.argv[1])
    lines = f.readlines()
    f.close()
 
    # File format:
    # [mismatch penalty]
    # [gap penalty]
    # [sequence 1]
    # [sequence 2]
 
    mismatch = int(lines[0])
    gap = int(lines[1])
    sequence1 = lines[2].strip()
    sequence2 = lines[3].strip()
 
    return mismatch, gap, sequence1, sequence2    
 
def needleman_wunsch(mismatch, gap, sequence1, sequence2):
    """
    Calculate the minimum penalty alignment.
    """
 
    # Get lengths of strings.
    n = len(sequence1)
    m = len(sequence2)
 
    # Make two-dimensional list for subproblem solutions.
    subproblems = [[0 for x in range(m+1)] for x in range(n+1)]
 
    # Fill in zeros on both dimensions with gap penalties.
    for i in range(n+1):
        subproblems[i][0] = i * gap
 
    for j in range(m+1):
        subproblems[0][j] = j * gap
 
    # Calculate subproblem solutions.
    for i in range(1, n+1):
        for j in range(1, m+1):
            case1 = subproblems[i-1][j-1]
 
            if sequence1[i-1] != sequence2[j-1]:
                case1 += mismatch
 
            case2 = subproblems[i-1][j] + gap
            case3 = subproblems[i][j-1] + gap
 
            subproblems[i][j] = min([case1, case2, case3])
 
    penalty = subproblems[n][m]
 
    # Backtrace to reconstruct optimal alignment.
    alignment1 = ""
    alignment2 = ""
 
    i = n
    j = m
    while i > 0 or j > 0:
        pos = subproblems[i][j]
        case1_match = subproblems[i-1][j-1]
        case1_mismatch = case1_match + mismatch
        case2 = subproblems[i-1][j] + gap
        case3 = subproblems[i][j-1] + gap
 
        if i > 0 and pos == case1_match:
            alignment1 = sequence1[i-1] + alignment1
            alignment2 = sequence2[j-1] + alignment2
            i -= 1
            j -= 1
        elif i > 0 and pos == case1_mismatch:
            alignment1 = sequence1[i-1] + alignment1
            alignment2 = sequence2[j-1] + alignment2
            i -= 1
            j -= 1
        elif i > 0 and pos == case2:
            alignment1 = sequence1[i-1] + alignment1
            alignment2 = '-' + alignment2
            i -= 1
        elif j > 0 and pos == case3:
            alignment1 = '-' + alignment1
            alignment2 = sequence2[j-1] + alignment2
            j -= 1
 
    return penalty, alignment1, alignment2
 
def main():
    """
    Read the penalties and sequences from the file,
    then calculate the optimal alignment.
    """
 
    mismatch, gap, sequence1, sequence2 = read()
 
    print "***Running Sequence Alignment***"
    print "Mismatch penalty:", mismatch
    print "Gap penalty:", gap
    print "Sequence 1 length:", len(sequence1)
    print "Sequence 2 length:", len(sequence2)
 
    penalty, alignment1, alignment2 = \
        needleman_wunsch(mismatch, gap, sequence1, sequence2)
 
    print "Total penalty:", penalty
    print "Optimal alignment"
    print
 
    line_length = 160
    for i in range(len(alignment1) / line_length + 1):
        start = i * line_length
        end = (i + 1) * line_length
        print alignment1[start:end]
        print alignment2[start:end]
        print
 
if __name__ == "__main__":
    main()