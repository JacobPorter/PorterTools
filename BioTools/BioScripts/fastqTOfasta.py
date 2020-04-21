#!/usr/bin/python
import optparse

from Bio import SeqIO


def main():

    p = optparse.OptionParser()

    p.add_option('--fastq', '-q', default='in.fastq', help='The fastq file.')
    p.add_option('--fasta', '-f', default='out.fa', help='The fasta file')
    options, _ = p.parse_args()

    fastq_loc = options.fastq
    fasta_loc = options.fasta

    input_handle = open(fastq_loc, "r")
    output_handle = open(fasta_loc, "w")

    sequences = SeqIO.parse(input_handle, "fastq")
    count = SeqIO.write(sequences, output_handle, "fasta")

    output_handle.close()
    input_handle.close()

    print("Converted fastq file %s to fasta file %s" % (fastq_loc, fasta_loc))
    print("Converted %i records" % count)


if __name__ == '__main__':
    main()
