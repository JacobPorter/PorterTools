#!/usr/bin/python
import SeqIterator
import datetime
import optparse
import os
import sys
import Constants
import re

logstr = "calculateSimulationAccuracy:\t"

def isUnmapped(flag):
    """Checks if the SAM flag indicates an unmapped read."""
    return ((int(flag) >> 2) % 2) == 1

def isFiltered(flag):
    """Checks if the filter bit is set."""
    return ((int(flag) >> 9) % 2) == 1

def isProper(flag):
    return ((int(flag) >> 1) % 2) == 1

def processSAM(SAM_file, gzip_switch, bound, paired_end):
    sam_iterator = SeqIterator.SeqIterator(SAM_file, file_type = Constants.SAM, gzip_switch = gzip_switch)
    counter_dict = {"Total_SAM_Records" : 0, "Records_Analyzed" : 0, "Correct" : 0, "Incorrect" : 0, "Unmap/Filter" : 0, "Improper" : 0} #"FP" : 0, "FN" : 0, "TP" : 0, "TN" : 0,
    record_dict = sam_iterator.convertToDict("R1", "R2")
    counter_dict["Total_SAM_Records"] = sam_iterator.records_processed()
    counter = 0
    for key in record_dict.keys():
        if counter < 10:
            sys.stderr.write("%s\n" % str(record_dict[key]))
            counter += 1
        if not paired_end and (len(record_dict[key]) > 1 or "XA:Z" in record_dict[key][0]):
#             if counter <= 10:
#                 print "0"
            continue
        if isUnmapped(record_dict[key][0][Constants.SAM_KEY_FLAG]) or isFiltered(record_dict[key][0][Constants.SAM_KEY_FLAG]):
#             if counter <= 10:
#                 print "1"
            counter_dict["Unmap/Filter"] += 1
            continue
        if paired_end and (len(record_dict[key])) != 2:
#             if counter <= 10:
#                 print "2"
            continue 
        if paired_end and ("XA:Z" in record_dict[key][0] or "XA:Z" in record_dict[key][1]):
#             if counter <= 10:
#                 print "3"
            continue
        if paired_end and (not isProper(record_dict[key][0][Constants.SAM_KEY_FLAG]) or not isProper(record_dict[key][1][Constants.SAM_KEY_FLAG])):
#             if counter <= 10:
#                 print "4"
            counter_dict["Improper"] += 1
            continue
        counter_dict["Records_Analyzed"] += 1
        SAM_record = record_dict[key][0]
        QNAME = re.split('_|:|-', SAM_record[Constants.SAM_KEY_QNAME])
        REAL_POS = int(QNAME[2])
        RNAME = SAM_record[Constants.SAM_KEY_RNAME]
        TEST_POS = int(SAM_record[Constants.SAM_KEY_POS])
        REAL_POS_u = REAL_POS + bound
        REAL_POS_l = REAL_POS - bound
        if paired_end:
            SAM_record2 = record_dict[key][1]
            RNAME2 = SAM_record2[Constants.SAM_KEY_RNAME]
            TEST_POS2 = int(SAM_record2[Constants.SAM_KEY_POS])
        if not paired_end and QNAME[1] == RNAME and  TEST_POS <= REAL_POS_u and TEST_POS >= REAL_POS_l:
            counter_dict["Correct"] += 1
        elif not paired_end:
            counter_dict["Incorrect"] += 1
        if paired_end and ((QNAME[1] == RNAME and TEST_POS <= REAL_POS_u and TEST_POS >= REAL_POS_l) or 
                           (QNAME[1] == RNAME2 and  TEST_POS2 <= REAL_POS_u and TEST_POS2 >= REAL_POS_l)):
            counter_dict["Correct"] += 1
        elif paired_end:
            counter_dict["Incorrect"] += 1
    return counter_dict


def report(counter_dict, paired_end, SAM_file, bounds, total_reads, now, later):
    total_reads = int(total_reads) + 0.0
    reads_analyzed = int(counter_dict["Records_Analyzed"]) + 0.0
    remaining_reads = total_reads - reads_analyzed
    add_s = "paired end" if paired_end else "single end" 
    if reads_analyzed == 0.0:
        reads_analyzed = 0.0000000000000001
    print "calculateSimulationAccuracy"
    print "---------------------------"
    print "For the %s SAM file %s, the process was started at %s and it took time %s." % (add_s, SAM_file, str(now), str(later - now))
    print "The interval around the real location had length %s above and below the real location." % str(bounds)
    print ""
    print "Total SAM records processed:\t%s" % (str(counter_dict["Total_SAM_Records"]))
    print "Total number of reads from the input:\t%s" % (str(total_reads))
    print "Total reads analyzed (uniquely mapped):\t%s" % (str(reads_analyzed))
    print "Reads unanalzyed, unmapped/filtered, improper, ambig:\t%s\t%s\t%s\t%s" % (str(remaining_reads), str(counter_dict["Unmap/Filter"]), str(counter_dict["Improper"]), str(remaining_reads - counter_dict["Unmap/Filter"] - counter_dict["Improper"]))
    print ""
    print "Percent correct of total    reads:\t%s" % (str(counter_dict["Correct"] / total_reads))
    print "Percent correct of analyzed reads:\t%s" % (str(counter_dict["Correct"] / reads_analyzed))
    print "Percent incorrect of total reads:\t%s" % (str(counter_dict["Incorrect"] / total_reads))
    print "Percent incorrect of reads analzyed:\t%s" % (str(counter_dict["Incorrect"] / reads_analyzed))
    print "Percent of reads         unanalyzed:\t%s" % (str(remaining_reads / total_reads))
    print ""
    print "Counter Object:"
    print str(counter_dict)


def main():
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <sam_file> <total_reads>"
    description = ""
    p = optparse.OptionParser(usage = usage, description = description)
    p.add_option('--gzip', '-z', help='The input file is gzip compressed, and the output file will be gzip compressed. [default: %default]', action='store_true', default=False)
    p.add_option('--bound', '-b', help='The interval length above and below the correct location in the reference genome.  This determines how sensitive the calculation is to the correct location. [default: %default]', default = 3)
    p.add_option('--paired_end', '-p', help='Turn this on if the data is paired end data.', action='store_true', default=False)
    options, args = p.parse_args()
    if len(args) == 0:
        p.print_help()
    if not os.path.exists(args[0]):
        p.error("The SAM file could not be found.")
    counter_dict = processSAM(args[0], options.gzip, int(options.bound), options.paired_end)
    later = datetime.datetime.now()
    report(counter_dict, options.paired_end, args[0], int(options.bound), int(args[1]), now, later)
    sys.stderr.write("%sThe process started at %s and took %s time.\n" % (logstr, str(now), str(later - now)))

if __name__ == "__main__":
    main()