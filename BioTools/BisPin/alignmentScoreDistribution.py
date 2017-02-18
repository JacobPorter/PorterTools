#!/usr/bin/python
import datetime
import optparse
import sys
import IndentedHelpFormatterWithNL 
import Constants
import SeqIterator
import math

bucket_count_default = 20

def findASDistribution(sam_file, bucket_count = bucket_count_default):
    sam_input = SeqIterator.SeqIterator(sam_file, file_type=Constants.SAM)
    as_max = -sys.maxint - 1
    as_min = sys.maxint
    no_hits = 0
    total_records = 0
    for record in sam_input:
        alignment_score = float(record[Constants.SAM_KEY_ALIGNMENT_SCORE])
        if record[Constants.SAM_KEY_RNAME].startswith(Constants.SAM_VALUE_STAR):
            no_hits += 1
            continue 
        if alignment_score > as_max:
            as_max = alignment_score
        if alignment_score < as_min:
            as_min = alignment_score
        total_records += 1
    #total_records = sam_input.records_processed()
    sam_input.reset()
    buckets = [0] * (bucket_count + 1)
    as_range = as_max - as_min
    bucket_length = as_range / bucket_count
    score_bounds = as_min
    score_bounds_list = []
    for _ in range(bucket_count):
        score_bounds += bucket_length
        score_bounds_list.append(score_bounds)
    stuff = 0
#     counter = 0
    total_score = 0
    median_position = int(total_records / 2)
    median_score = 0
    for record in sam_input:
        alignment_score = float(record[Constants.SAM_KEY_ALIGNMENT_SCORE])
        if record[Constants.SAM_KEY_RNAME].startswith(Constants.SAM_VALUE_STAR):
            continue
        if stuff == median_position:
            median_score = alignment_score
        total_score += alignment_score
        bucket_index = int(math.floor((alignment_score - as_min) / bucket_length))
#         if bucket_index == 2:
#             print bucket_index, alignment_score
        buckets[bucket_index] += 1
#         if counter < 10:
#             print bucket_index, alignment_score
#             counter += 1
        stuff += 1
    avg_score = total_score / (stuff + 0.0)
    assert total_records == stuff
    freq = map((lambda x: x / (total_records + 0.0)), buckets)
    return (freq, score_bounds_list, total_records, as_max, as_min, buckets, avg_score, median_score, no_hits)

def main():
    """
    Calculates a distribution of the alignment scores for a SAM file.  The distribution is printed to standard out.
    """
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <sam_file> "
    version = "%prog " + Constants.version
    description = ""
    epilog = ""
    p = optparse.OptionParser(usage = usage, version = version, 
                              description = description, epilog = epilog, 
                              formatter = IndentedHelpFormatterWithNL.IndentedHelpFormatterWithNL())
    p.add_option('--buckets','-b',help='The number of buckets to produce in the frequency distribution. [default: %default]', default = bucket_count_default)
    options, args = p.parse_args()
    if len(args) != 1:
        p.print_help()
        p.error("There must be one argument.")
    sam_file = args[0]
    sys.stdout.write("Processing the file %s with %s buckets.\n" % (sam_file, str(options.buckets)))
    stats = findASDistribution(sam_file, bucket_count = int(options.buckets))
    (freq, score_bounds_list, total_records, as_max, as_min, buckets, avg_score, median_score, no_hits) = stats
    stats = map(lambda x: str(x), stats)
    #print len(stats)
    #sys.stdout.write(str(stats[0]) + "\n" + str(stats[1]) + "\n")
    result_string = "Frequency distribution: %s\nScore bounds: %s\nReads processed: %s\nMaximum score: %s\nMinimum score: %s\nCounts: %s\nAverage score: %s\nMedian score: %s\nNumber of no hit alignments: %s\n" % (freq, score_bounds_list, total_records, as_max, as_min, buckets, avg_score, median_score, no_hits)
    later = datetime.datetime.now()
    runtime = later - now
    sys.stdout.write(result_string + "Runtime: " + str(runtime) + "\n")

if __name__ == "__main__":
    main()