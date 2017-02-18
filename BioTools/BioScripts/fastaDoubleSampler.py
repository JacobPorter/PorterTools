#!/usr/bin/python
import optparse, math, random, sys, itertools
import SeqGenerator
from os import listdir
from os.path import isfile, join

'''

'''

def main():
    p = optparse.OptionParser()
    p.add_option('--fasta_dir', '-d', help='A directory of fasta files to sample from')
    p.add_option('--fasta_dir2', '-t', help='A directory of fasta files that match up with the other fasta file directory.  The files must be sorted to work.  The files must be sortable in the same way.')
    p.add_option('--output', '-o', help='The prefix for the output files.  Each file will have a number appended to it indicating the sample number')
    p.add_option('--output2','-u',help='The prefix for the output files for the second set of fasta files.')
    p.add_option('--percent', '-p', default='10.0', help='Percent of the total sample to take.  Must be a real value between 0 and 100')
    p.add_option('--num', '-n', default='10', help='Number of samples to take.')
    options, _ = p.parse_args()
    
    sampler(options.fasta_dir, options.fasta_dir2, options.output, options.output2, options.percent, options.num)
    print "Done!"
    
def sampler(fasta_dir, fasta_dir2, output, output2, percent, num):
    num = int(num)
    percent = float(percent)
    if percent < 0 or percent > 100:
        raise ValueError
    percent = percent / 100.0
    onlyfiles = [ f for f in listdir(fasta_dir) if isfile(join(fasta_dir,f)) ]
    onlyfiles.sort()
    
    onlyfiles2 = [ f for f in listdir(fasta_dir2) if isfile(join(fasta_dir2,f)) ]
    onlyfiles2.sort()
    #onlyfiles = tuple(onlyfiles)
    print "Counting the sequences..."
    sys.stdout.flush()
    counts = count_sequences(fasta_dir, onlyfiles)
    total_seqs = sum(counts)
    num_to_get = int(math.ceil(percent * total_seqs))
    print "From a total of %s sequences in %s files, %s sequences will be sampled uniformly at random with replacement for each file." % (str(total_seqs) , str(len(onlyfiles)), str(num_to_get))
    sys.stdout.flush()
    locations = []
    print "Generating random numbers..."
    sys.stdout.flush()
    for i in range(num):
        for j in range(int(num_to_get)):
            r_n = random.randint(1, total_seqs)
            s = (r_n, i)
            locations.append(s)
    print "Done generating random numbers."
    sys.stdout.flush()
    locations.sort(key=lambda tup: tup[0])
    files = []
    files2 = []
    print "Creating %s files with prefix %s with a total of %s random sequences with replacement to be drawn." % (str(num), output, str(len(locations)))
    for i in range(num):
        files.append(SeqGenerator.SeqWriter(open(join(output + "." + str(i) + ".fa"), 'w' )))
        files2.append(SeqGenerator.SeqWriter(open(join(output2 + "." + str(i) + ".fa"), 'w' )))
    seq_counter = 0
    list_counter = 0
    for ii in range(len(onlyfiles)):
        f = join(fasta_dir, onlyfiles[ii])
        f2 = join(fasta_dir2, onlyfiles2[ii])
        seq_g = SeqGenerator.SeqGenerator(f)
        seq_g2 = SeqGenerator.SeqGenerator(f2)
        print "Processing %s and %s with %s sequences" % (f, f2, str(counts[ii]))
        sys.stdout.flush()
        for seq,seq2 in itertools.izip(seq_g, seq_g2):
            seq_counter += 1
            while list_counter < len(locations) and seq_counter == locations[list_counter][0]:
                files[locations[list_counter][1]].write(seq)
                files2[locations[list_counter][1]].write(seq2)
                list_counter += 1
            if list_counter >= len(locations):
                return
            if list_counter % 10000 == 0:
                for i in range(num):
                    files[i].flush()
                    files2[i].flush()
    
    
    
def count_sequences(my_dir, onlyfiles):
    count = []
    for f in onlyfiles:
        print f
        sys.stdout.flush()
        seq_g = SeqGenerator.SeqGenerator(join(my_dir, f))
        count.append(seq_g.count())
    return tuple(count)

if __name__ == '__main__':
    main()