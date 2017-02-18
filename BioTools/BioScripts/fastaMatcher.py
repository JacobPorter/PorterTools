#!/usr/bin/python
import optparse, sys
import SeqGenerator
from os import listdir
from os.path import isfile, join

'''
Scan through all the sequences in the files in one directory and then match them up with all the sequences in the files in another directory.  The sequence ids must match.
'''

def main():
    p = optparse.OptionParser()
    p.add_option('--match_to_dir', '-m', help='A directory of fasta files to get matches for')
    p.add_option('--match_from_dir', '-f', help='A directory of fasta files to get matches from.  Assumes the sequences are sorted.')
    p.add_option('--output', '-o', help='The prefix for the output files.  Each file will have a number appended to it indicating the sample number')
    options, _ = p.parse_args()
    
    matcher(options.match_to_dir, options.match_from_dir, options.output)
    print "Done!"
    
def matcher(match_to_dir, match_from_dir, output):
    onlyFromFiles = [ f for f in listdir(match_from_dir) if isfile(join(match_from_dir,f)) ]
    onlyFromFiles.sort()
    onlyToFiles = [ f for f in listdir(match_to_dir) if isfile(join(match_to_dir,f)) ]
    onlyToFiles.sort()
    seq_ids = get_sequence_ids(match_to_dir, onlyToFiles)
    seq_ids.sort(key=lambda tup: tup[0])
    seq_ids = iter(seq_ids)
    
    locations = {}
    for f in onlyToFiles:
        locations[f] = SeqGenerator.SeqWriter(open(join(output + "." + f + ".fa"), 'w' ))
    
    count = 0
    for f in onlyFromFiles:
        try:
            seq_to = seq_ids.next()
        except StopIteration:
            break
        seq_g = SeqGenerator.SeqGenerator(join(match_from_dir, f))
        for s in seq_g:
            if str(s[0].split('#')[0]) == str(seq_to[0].split('#')[0]):
                locations[seq_ids[count][1]].write(s)
                count += 1
                try:
                    seq_to = seq_ids.next()
                except StopIteration:
                    break
            if count >= 1000000:
                count = 0
                for f1 in onlyToFiles:
                    locations[f1].flush()
    
    

def get_sequence_ids(my_dir, onlyfiles):
    seq_ids = []
    for f in onlyfiles:
        print f
        sys.stdout.flush()
        seq_g = SeqGenerator.SeqGenerator(join(my_dir, f))
        for seq in seq_g:
            seq_ids.append((seq[0], f))
    return seq_ids

if __name__ == '__main__':
    main()