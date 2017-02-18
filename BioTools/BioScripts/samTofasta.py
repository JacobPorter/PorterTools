#!/usr/bin/python
import sys, optparse

def main():
    p = optparse.OptionParser()
    
    p.add_option('--sam_file', '-s' , help='The sam input file.')
    p.add_option('--fasta_file', '-f', help='The fasta file to be written to.  If the file exists, it will be overwritten.')
    p.add_option('--ignore_string', '-i', help='Ignores sam entries with the user specified string in the second column.  For example, bsmap puts an asterisk in the second column for some sequences.')
    options, _ = p.parse_args()
    
    convertSamTOFasta(options.sam_file, options.fasta_file, options.ignore_string)


def write_fasta(fasta_fd, lst):
    fasta_fd.write('>' + lst[0] + '\n')
    fasta_fd.flush()
    fasta_fd.write(lst[9] + '\n')
    fasta_fd.flush()
    

def convertSamTOFasta(sam_fn, fasta_fn, ignore_string):
    sam_fd = open(sam_fn)
    fasta_fd = open(fasta_fn, 'w')
    count = 0
    for line in sam_fd:
        if line[0] != '@':
            my_spl = line.split('\t')
            if ignore_string == None:
                write_fasta(fasta_fd, my_spl)
                count += 1
            elif my_spl[2].strip() != ignore_string:
                write_fasta(fasta_fd, my_spl)
                count += 1
    sys.stdout.write("Created " + str(count) + " sequences from sam file " + sam_fn + " stored to the fasta file " + fasta_fn + "\n")
    sys.stdout.flush()
    sam_fd.close()
    fasta_fd.close()
            
            
if __name__ == '__main__':
    main()