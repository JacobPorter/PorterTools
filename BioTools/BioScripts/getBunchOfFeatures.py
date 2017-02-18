#!/usr/bin/python
import optparse
import sys
import os

def main():
    p = optparse.OptionParser(usage="usage: %prog [options] path_to_feature_files", 
                              version="%prog 0.0")
    _, args = p.parse_args()
    if len(sys.argv[1:]) == 0:
        print "%sA required argument is missing.  Use -h or --help for argument information."
        p.print_help()
        exit(1)
    mypath = args[0]
    onlyfiles = [ f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath,f)) ]
    onlyfiles.sort()
    matrix = []
    get_classes = True
    classes = []
    for f in onlyfiles:
        count = 0
        feature_vector = []
        for line in open(os.path.join(mypath, f)):
            if count <= 1:
                pass
            else:
                try:
                    feature_vector.append(float(line.split("\t")[0]))
                    if get_classes:
                        classes.append(int(line.strip("\n").split("\t")[1]))
                except ValueError:
                    pass
            count += 1
        get_classes = False
        matrix.append(feature_vector)
    print str(len(matrix[0]))
    print str(len(matrix))
    for i in xrange(len((matrix[0]))):
        str1 = ""
        for j in xrange(len(matrix)):
            str1 += str(matrix[j][i]) + "\t"
        str1 += str(classes[i])
        print str1
        
    
if __name__ == '__main__':
    main()