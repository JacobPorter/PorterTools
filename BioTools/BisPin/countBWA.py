#!/usr/bin/python
import SeqIterator
import datetime
import optparse
import os
import sys
import Constants

logstr = "countBWA:\t"

"""
A program that counts how many ambiguously mapped, unmapped, and uniquely mapped reads
there are from a BWA produced SAM file.
@author: Jacob Porter
"""

class Counts:
    """
    A class that stores the counts.
    """
    def __init__(self):
        self.__total = 0
        self.__ambig = 0
        self.__unmap = 0
        self.__filt = 0
        self.__uniq = 0
        self.__sam = 0
    
    def count_sam(self):
        self.__sam += 1
        
    def count_total(self):
        self.__total += 1
    
    def count_ambig(self):
        self.__ambig += 1
        
    def count_unmap(self):
        self.__unmap += 1
        
    def count_uniq(self):
        self.__uniq += 1
        
    def count_filt(self):
        self.__filt += 1    
    
    def total(self):
        return self.__total
    
    def ambig(self):
        return self.__ambig
    
    def unmap(self):
        return self.__unmap
    
    def filt(self):
        return self.__filt
    
    def uniq(self):
        return self.__uniq
    
    def sam(self):
        return self.__sam
    

def isUnmapped(flag):
    """Checks if the SAM flag indicates an unmapped read."""
    return ((int(flag) >> 2) % 2) == 1

def isFiltered(flag):
    """Checks if the filter bit is set."""
    return ((int(flag) >> 9) % 2) == 1

def isProper(flag):
    return ((int(flag) >> 1) % 2) == 1

def processSAM(SAM_file, gzip_switch, paired_end):
    sam_iterator = SeqIterator.SeqIterator(SAM_file, file_type = Constants.SAM, gzip_switch = gzip_switch)
    my_counter = Counts()
    record_dict = sam_iterator.convertToDict("R1", "R2")
    my_counter.__sam = sam_iterator.records_processed()
    if not paired_end:
        processSingle(record_dict, my_counter)
    else:
        processPaired(record_dict, my_counter)
    return my_counter

def processSingle(record_dict, my_counter):
    for key in record_dict:
        SAM_record_list = record_dict[key]
        my_counter.count_total()
        if len(SAM_record_list) != 1:
            my_counter.count_ambig()
        elif None != SAM_record_list[0].get("XA:Z", None):
            my_counter.count_ambig()
        elif isUnmapped(SAM_record_list[0][Constants.SAM_KEY_FLAG]):
            my_counter.count_unmap()
        elif isFiltered(SAM_record_list[0][Constants.SAM_KEY_FLAG]):
            my_counter.count_filt()
        else:
            my_counter.count_uniq()

def processPaired(record_dict, my_counter):
    for key in record_dict:
        SAM_record_list = record_dict[key]
        my_counter.count_total()
        if len(SAM_record_list) != 2:
            my_counter.count_ambig()
        elif (None != SAM_record_list[0].get("XA:Z", None)) or (None != SAM_record_list[1].get("XA:Z", None)):
            my_counter.count_ambig()
        elif isUnmapped(SAM_record_list[0][Constants.SAM_KEY_FLAG]) or isUnmapped(SAM_record_list[1][Constants.SAM_KEY_FLAG]):
            my_counter.count_unmap()
        elif not isProper(SAM_record_list[0][Constants.SAM_KEY_FLAG]) or not isProper(SAM_record_list[1][Constants.SAM_KEY_FLAG]):
            my_counter.count_unmap()
        elif isFiltered(SAM_record_list[0][Constants.SAM_KEY_FLAG]) or isFiltered(SAM_record_list[1][Constants.SAM_KEY_FLAG]):
            my_counter.count_filt()
        else:
            my_counter.count_uniq()

def report(counter, paired_end, SAM_file, now, later):
    add_s = "paired end" if paired_end else "single end" 
    sys.stdout.write("CountBWA by Jacob Porter\n")
    sys.stdout.write("\n")
    sys.stdout.write("The %s SAM file %s had %s SAM records.  Processing began at %s taking %s time.\n" % (add_s, str(SAM_file), str(counter.sam()), str(now), str(later - now)))
    t = counter.total() + 0.0
    if t == 0.0:
        t = 0.000000000001
    sys.stdout.write("Reads:\t%s\n" % str(counter.total()))
    sys.stdout.write("Unique:\t%s\t%f\n" % (str(counter.uniq()), counter.uniq() / t))
    sys.stdout.write("Ambig:\t%s\t%f\n" % (str(counter.ambig()), counter.ambig() / t))
    sys.stdout.write("Unmap:\t%s\t%f\n" % (str(counter.unmap()), counter.unmap() / t))
    sys.stdout.write("Filt:\t%s\t%f\n" % (str(counter.filt()), counter.filt() / t))
    sys.stdout.write("Unmap + Filt:\t%s\t%f\n" % (str(counter.filt() + counter.unmap()), (counter.filt() + counter.unmap()) / t))
    sys.stdout.flush()

def main():
    now = datetime.datetime.now()
    usage = "usage: %prog [options] <sam_file> "
    description = ""
    p = optparse.OptionParser(usage = usage, description = description)
    p.add_option('--gzip', '-z', help='The input file is gzip compressed, and the output file will be gzip compressed. [default: %default]', action='store_true', default=False)
    p.add_option('--paired_end', '-p', help='Turn this on if the data is paired end data.', action='store_true', default=False)
    p.add_option('--bound', '-b', help='The interval length above and below the correct location in the reference genome.  This determines how sensitive the calculation is to the correct location. [default: %default]', default = 3)
    options, args = p.parse_args()
    if len(args) == 0:
        p.print_help()
    if not os.path.exists(args[0]):
        p.error("The SAM file could not be found.")
    counter = processSAM(args[0], options.gzip, options.paired_end)
    later = datetime.datetime.now()
    report(counter, options.paired_end, args[0], now, later)
    sys.stderr.write("%sThe process started at %s and took %s time.\n" % (logstr, str(now), str(later - now)))

if __name__ == "__main__":
    main()