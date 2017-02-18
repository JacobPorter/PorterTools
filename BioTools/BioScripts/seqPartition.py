#!/usr/bin/python
import SeqGenerator
import optparse, sys, math, os


def main():
    p = optparse.OptionParser()
    p.add_option('--input','-i',help='Sequence file to partition')
    p.add_option('--output', '-o', help='File prefix for writing to.  File will be overwritten if it exists.')
    p.add_option('--percent','-p',help='Percent of sequences in each file.  Must be a number between 0 and 100')
    p.add_option('--type', '-t', help='The file type to partition.  Must be fasta or fastq.', default='fasta')
    options,_ = p.parse_args()
    seq_part(options.input, options.output, options.type, options.percent)
    print "Done!"
    
def seq_part(file_in, file_out, file_type, percent):
    percent = float(percent)
    if percent < 0 or percent > 100:
        print "Percent must be a number between 0 and 100"
        raise ValueError
    percent = percent / 100.0
    print "Sequence file to partition %s of type %s" % (file_in, file_type)
    sys.stdout.flush()
    seq_g = SeqGenerator.SeqGenerator(file_in, file_type=file_type)
    count = seq_g.count()
    my_number = int(math.ceil(count * percent))
    print "Each partition will have about %d sequences, which is %s percent" % (my_number, str(percent))
    sys.stdout.flush()
    seq_g = SeqGenerator.SeqGenerator(file_in, file_type=file_type)
    part_cnt = 0
    part = 0
    seq_w = SeqGenerator.SeqWriter(open(os.path.join(file_out + "." + str(part) + "." + file_type), 'w'), file_type=file_type)
    for s in seq_g:
#         if part_cnt == 0:
#             print s
        if part_cnt == my_number:
            seq_w.flush()
            part += 1
            part_cnt = 0
            seq_w = SeqGenerator.SeqWriter(open(os.path.join(file_out + "." + str(part) + "." + file_type), 'w'), file_type=file_type)
        seq_w.write(s)
        part_cnt += 1
    seq_w.flush()        
            
    
if __name__ == '__main__':
    main()