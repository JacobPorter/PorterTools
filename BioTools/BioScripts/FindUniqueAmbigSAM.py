#!/usr/bin/python
import optparse, sys

def main():
    p = optparse.OptionParser()
    
    p.add_option('--file', '-f' , help='The SAM file to calculate ambig, etc from.')
    p.add_option('-v', action='store_true', dest='verbose', help='Prints verbose output')
    options, _ = p.parse_args()
    
    print "Running FindUniqueAmbigSAM.py on SAM file %s" % options.file
    unique_c = 0
    total_c = 0
    ambig_c = 0
    unmap_c = 0
    sam_file = open(options.file)
    sys.stdout.flush()
    for line in sam_file:
        if line[0] != '@':
            total_c += 1
            entry_lst = line.split()
            AS = entry_lst[11]
            if AS[0:2] != 'AS':
                if options.verbose:
                    print "ERROR!  AS not found: \n %s" % (line)
                    sys.stdout.flush()
                unmap_c += 1
                continue
            XS = entry_lst[12]
            AS_n = AS.split('i:')[1]
            XS_n = XS.split('i:')[1]
            if XS[0:2] != 'XS' or AS_n != XS_n:
                unique_c += 1
                if options.verbose and XS[0:2] != 'XS':
                    print "ERROR!  XS not found: \n %s" % (line)
                    sys.stdout.flush()
            elif AS_n == XS_n:
                ambig_c += 1
        else:
            continue
    print "Total reads aligned: %d" % (total_c)
    print "Unique reads: %d" % (unique_c)
    print "Ambiguous reads: %d" % (ambig_c) 
    print "Unmapped reads: %d" % (unmap_c)
    print "Done!"

if __name__ == '__main__':
    main()