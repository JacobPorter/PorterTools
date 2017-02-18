#!/usr/bin/python
import optparse
import sys
import os

def main():
    p = optparse.OptionParser(usage="usage: %prog [options] path_to_outputs", 
                              version="%prog 0.0")
    _, args = p.parse_args()
    if len(sys.argv[1:]) == 0:
        print "%sA required argument is missing.  Use -h or --help for argument information."
        p.print_help()
        exit(1)
    mypath = args[0]
    onlyfiles = [ f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath,f)) ]
    onlyfiles.sort()
    performance_dict = []
    for f in onlyfiles:
        performance_list = [f]
        for line in open(os.path.join(mypath, f)):
            if line.startswith("Recognition rate:"):
                performance_list.append(float(line.strip("\n").split(":")[1].replace("%", "")))
            elif line.startswith("C1:"):
                performance_list.append(float(line.strip("\n").split(":")[1]))
            elif line.startswith("C2:"):
                performance_list.append(float(line.strip("\n").split(":")[1]))
            elif line.startswith("C3:"):
                performance_list.append(float(line.strip("\n").split(":")[1]))
        performance_dict.append(performance_list)
    performance_dict.sort(key = lambda row: row[1])
    print "Feature \t Recognition rate \t C1 \t C2 \t C3"
    for row in performance_dict:
        str1 = ""
        for item in row:
            str1 += str(item) + "\t"
        print str1
    
if __name__ == '__main__':
    main()
