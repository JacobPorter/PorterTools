#!/usr/bin/python
import SeqGenerator
import optparse, sys
from os import listdir, remove, rename
from os.path import isfile, join

def main():
    p = optparse.OptionParser()
    p.add_option('--input','-i',help='Directory of fasta files to uniquify.  Destructive.')
    options,_ = p.parse_args()
#    if len(args) == 0:
#        print "Use --help for help documentation"
#        return
    print "Running fastaUniquify on fasta file directory: %s" % options.input
    sys.stdout.flush()
    fastaUnique(options.input)
    print "Done!"

def fastaUnique(fasta_dir):
    onlyfiles = [ f for f in listdir(fasta_dir) if isfile(join(fasta_dir,f)) ]
    onlyfiles.sort()
    for f in onlyfiles:
        print "Doing file %s" % f
        sys.stdout.flush()
        temp_fd = open(join(fasta_dir, "temp.fq"), 'w')
        seq_w = SeqGenerator.SeqWriter(temp_fd, file_type="fastq")
        seq_r = SeqGenerator.SeqGenerator(join(fasta_dir, f), file_type="fastq")
        last_id = ""
        count = 1
        for s in seq_r:
            parts = s[0].split(' ')
            this_id = parts[0]
            the_rest = ""
            for i in range(len(parts)):
                if i == 0:
                    continue
                the_rest += " " + parts[i]
            if this_id == last_id:
                temp_id = this_id + "_" +  str(count) + the_rest
                seq_w.write((temp_id, s[1], s[2]))
                count += 1
            else:
                count = 1
                last_id = this_id
                seq_w.write(s)
        seq_w.flush()
        seq_w = None
        temp_fd.close()
        remove(join(fasta_dir, f))
        rename(join(fasta_dir, "temp.fq"), join(fasta_dir, f))

if __name__ == '__main__':
    main()