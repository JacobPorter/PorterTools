#!/usr/bin/python
from Bio import SeqIO
import optparse

def main():

    p = optparse.OptionParser()

    p.add_option('--fastq', '-q', default='out.fastq', help='The fastq file to convert to a fasta file.')
    p.add_option('--fasta', '-f', default='in.fa', help='The fasta file to save the output.  Will be overwritten.')
    options, _ = p.parse_args()

    fastq_loc= options.fastq
    fasta_loc = options.fasta


    input_handle = open(fasta_loc, "r")
    output_handle = open(fastq_loc, "w")

    sequences = list(SeqIO.parse(input_handle, "fasta"))
    for seq in sequences:
        seq.letter_annotations["phred_quality"] = [26] * len(seq)
    count = SeqIO.write(sequences, output_handle, "fastq")

    output_handle.close()
    input_handle.close()

    print("Converted fastq file %s to fasta file %s" % (fastq_loc, fasta_loc))
    print("Converted %i records" % count)

if __name__ == '__main__':
    main()
