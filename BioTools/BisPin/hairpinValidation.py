#!/usr/bin/python
import datetime
import optparse
import os
import sys
import SeqIterator
import Constants

logstr = "hairpinValidation:\t"

def isUnmapped(flag):
    """Checks if the SAM flag indicates an unmapped read."""
    return ((int(flag) >> 2) % 2) == 1

def isFiltered(flag):
    """Checks if the filter bit is set."""
    return ((int(flag) >> 9) % 2) == 1

def hairpinValidate(hairpin_recovered_filename, test_filename, bound, gzip_switch):
    counter = {"Hairpin_Hits" : 0, "Test_NotFound" : 0, "Test_Ambig" : 0, "Test_UnmapFilt" : 0, "Test_Correct" : 0, "Test_Incorrect" : 0}
    hairpin_dict = SeqIterator.SeqIterator(hairpin_recovered_filename, file_type = Constants.SAM).convertToDict()
    test_dict = SeqIterator.SeqIterator(test_filename, file_type = Constants.SAM, gzip_switch = gzip_switch).convertToDict()
    for H_QNAME in hairpin_dict:
        if len(hairpin_dict[H_QNAME]) != 1 or "XA:Z" in hairpin_dict[H_QNAME][0]:
            continue
        hairpin_SAM = hairpin_dict[H_QNAME][0]
        flag = hairpin_SAM[Constants.SAM_KEY_FLAG]
        if isUnmapped(flag) or isFiltered(flag):
            continue
        counter["Hairpin_Hits"] += 1
        hairpin_pos = int(hairpin_SAM[Constants.SAM_KEY_POS])
        test_SAM_list = test_dict.get(H_QNAME, None)
        if test_SAM_list == None:
            counter["Test_NotFound"] += 1
        elif len(test_SAM_list) != 1 or "XA:Z" in test_SAM_list[0]:
            counter["Test_Ambig"] += 1
        else:
            test_pos = int(test_SAM_list[0][Constants.SAM_KEY_POS])
            test_SAM = test_SAM_list[0]
            if isUnmapped(test_SAM[Constants.SAM_KEY_FLAG]) or isFiltered(test_SAM[Constants.SAM_KEY_FLAG]):
                counter["Test_UnmapFilt"] += 1
            elif test_SAM[Constants.SAM_KEY_RNAME] == hairpin_SAM[Constants.SAM_KEY_RNAME] and \
                test_pos <= hairpin_pos + bound and \
                test_pos >= hairpin_pos - bound:
                counter["Test_Correct"] += 1
            else:
                counter["Test_Incorrect"] += 1
    return counter

def report(counter, hairpin_recovered_filename, test_filename, bound, now, later):
    reads_analyzed = int(counter["Hairpin_Hits"]) + 0.0
    if reads_analyzed == 0.0:
        reads_analyzed = 0.0000000000000001
    print "hairpinValidation"
    print "---------------------------"
    print "Using the %s hairpin recovery SAM file, processing the SAM file %s took time %s beginning at %s." % (hairpin_recovered_filename, test_filename, str(later - now), str(now))
    print "The interval around the real location had length %s above and below the real location." % str(bound)
    print ""
    keys = counter.keys()
    keys.sort()
    for key in keys:
        print "%s reads:\t%d\t%f" % (key, counter[key], counter[key] / reads_analyzed)

def main():
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <hairpin_recovered_sam_file> <test_sam_file> "
    description = ""
    p = optparse.OptionParser(usage = usage, description = description)
    p.add_option('--gzip', '-z', help='The input test SAM file is gzip compressed, and the output file will be gzip compressed. [default: %default]', action='store_true', default=False)
    #p.add_option('--paired_end', '-p', help='Turn this on if the data is paired end data.', action='store_true', default=False)
    p.add_option('--bound', '-b', help='The interval length above and below the correct location in the reference genome.  This determines how sensitive the calculation is to the correct location. [default: %default]', default = 3)
    options, args = p.parse_args()
    if len(args) == 0:
        p.print_help()
    if len(args) != 2:
        p.error("There must be two arguments given.")
    if not os.path.exists(args[0]):
        p.error("The hairpin recovered SAM file could not be found.")
    if not os.path.exists(args[1]):
        p.error("The test SAM file could not be found.")
    try:
        bound = int(options.bound)
        if bound < 0:
            raise ValueError
    except ValueError:
        p.error("The bounds must be an integer larger or equal to zero.")
    counter = hairpinValidate(args[0], args[1], bound, options.gzip)
    later = datetime.datetime.now()
    report(counter, args[0], args[1], bound, now, later)
    sys.stderr.write("%sThe process started at %s and took %s time.\n" % (logstr, str(now), str(later - now)))

if __name__ == "__main__":
    main()