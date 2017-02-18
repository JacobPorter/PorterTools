#!/usr/bin/python
import optparse, math, random, sys
import SeqIterator
from os import listdir
from os.path import isfile, join

'''
Utility program for generating random replicates with replacement from a directory of sequence files.
'''

def main():
    p = optparse.OptionParser()
    p.add_option('--dir', '-d', help='A directory of files to sample from.  There must only be files of the appropriate type in the directory.')
    p.add_option('--output', '-o', help='The prefix for the output files.  Each file will have a number appended to it indicating the sample number')
    p.add_option('--num', '-n', default='.10', help='Total number of records to sample for each replicate.  If the value is larger than 0, then that number of records will be sampled from each replicate.  If the value is between 0 and 1 non inclusive, then that percent of the total number of records will be sampled for each replicate.' )
    p.add_option('--replicates', '-r', default='10', help='Number of replicates to generate.')
    p.add_option('--type', '-t', default='fastq', help="The type of files to sample from.  Must be either 'fasta', 'fastq', or 'sam'")
    options, _ = p.parse_args()
    sampler(options.dir, options.output, float(options.num), int(options.replicates), options.type)
    print "Done!"
    
    
def sampler(the_dir, output, num, replicates, file_type):
    if num < 0:
        raise ValueError
    onlyfiles = [ f for f in listdir(the_dir) if isfile(join(the_dir,f)) ]
    onlyfiles.sort()
    #onlyfiles = tuple(onlyfiles)
    print "Counting the sequences..."
    sys.stdout.flush()
    counts = count_sequences(the_dir, onlyfiles, file_type)
    total_seqs = sum(counts)
    if num < 1:
        num_to_get = int(math.ceil(num * total_seqs))
    else:
        num_to_get = num
    print "From a total of %s sequences in %s files, %s sequences will be sampled uniformly at random with replacement for each of %s files." % (str(total_seqs) , str(len(onlyfiles)), str(num_to_get), str(replicates))
    sys.stdout.flush()
    locations = []
    for i in range(replicates):
        for _ in range(int(num_to_get)):
            r_n = random.randint(1, total_seqs)
            s = (r_n, i)
            locations.append(s)
    locations.sort(key=lambda tup: tup[0])
    files = []
    print "Creating %s files with prefix %s with a total of %s random sequences with replacement to be drawn." % (str(replicates), output, str(len(locations)))
    for i in range(replicates):
        files.append(SeqIterator.SeqWriter(open(join(output + "." + str(i) + "." + file_type), 'w' ), file_type=file_type ))
    seq_counter = 0
    list_counter = 0
    for ii in range(len(onlyfiles)):
        f = join(the_dir, onlyfiles[ii])
        seq_g = SeqIterator.SeqIterator(f, file_type=file_type)
        print "Processing %s with %s sequences" % (f, str(counts[ii]))
        sys.stdout.flush()
        for seq in seq_g:
            seq_counter += 1
            while list_counter < len(locations) and seq_counter == locations[list_counter][0]:
                files[locations[list_counter][1]].write(seq)
                list_counter += 1
            if list_counter >= len(locations):
                return
            if list_counter % 10000 == 0:
                for i in range(replicates):
                    files[i].flush()
    
    
def count_sequences(my_dir, onlyfiles, file_type):
    count = []
    for f in onlyfiles:
        print f
        sys.stdout.flush()
        seq_g = SeqIterator.SeqIterator(join(my_dir, f), file_type)
        count.append(seq_g.count())
    return tuple(count)

if __name__ == '__main__':
    main()