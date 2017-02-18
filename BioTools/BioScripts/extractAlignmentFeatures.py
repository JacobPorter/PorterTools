#!/usr/bin/python
import optparse
import sys
import SeqGenerator
import csv
import math
import nw_hairpin

def main():
    p = optparse.OptionParser(usage="usage: %prog [options] CtoT_file GtoA_file sam_file ambig_file unmap_file feature_file ", 
                              version="%prog 0.0")
    options, args = p.parse_args()
    if len(sys.argv[1:]) == 0:
        print "%sA required argument is missing.  Use -h or --help for argument information."
        p.print_help()
        exit(1)    
    CtoT_file = args[0]
    GtoA_file = args[1]
    path_to_feature_file = args[5]
    extractFeatures(CtoT_file, GtoA_file, path_to_feature_file, args[2], args[3], args[4])
    
    
def extractFeatures(CtoT_file, GtoA_file, path_to_feature_file, sam_file, ambig_file, unmap_file):
    sam_seq_g = SeqGenerator.SeqGenerator(sam_file, file_type="SAM")
    align_type = {}
    for record in sam_seq_g:
        align_type[record["seq_id"]] = 1 #"unique"
    ambig_seq_g = SeqGenerator.SeqGenerator(ambig_file)
    for record in ambig_seq_g:
        align_type[record[0]] = 2#"ambig"
    unmap_seq_g = SeqGenerator.SeqGenerator(unmap_file)
    for record in unmap_seq_g:
        align_type[record[0]] = 3#"unmap"
    GtoA_gen = SeqGenerator.SeqGenerator(GtoA_file)
    GtoA_sequences = {}
    for record in GtoA_gen:
        GtoA_sequences[record[0].split('#')[0]] = record[1]
    feature_file_desc = open(path_to_feature_file, 'w')
    feature_writer = csv.writer(feature_file_desc, delimiter='\t',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    feature_header = ["length", "seq_entropy", "first_20_entropy", 
                    "A_freq", "C_freq",  "G_freq", "T_freq", "N_freq",
                    "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                    "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", 
                    "length_X_seq_entropy", "mismatch_gap_percent" , "mapping_type"]
    #feature_writer.writerow(feature_header)
    seq_g = SeqGenerator.SeqGenerator(CtoT_file)
    for record in seq_g:
        DNA_seq = record[1].upper()
        feature_row = []
        #feature_row = [record[0]]
        feature_row.append(len(DNA_seq))
        DNA_entropy = findEntropy(DNA_seq)
        DNA_first_20_entropy = findEntropy(DNA_seq[0:20])
        feature_row.append(DNA_entropy[0])
        feature_row.append(DNA_first_20_entropy[0])
        for f in DNA_entropy[1]:
            feature_row.append(f)
        dimer_counts = []
        for dimer in feature_header[9:26]:  # Check features
            dimer_counts.append(countKmers(DNA_seq, dimer))
        total_dimers = sum(dimer_counts)
        for dimer_c in dimer_counts:
            feature_row.append(dimer_c / (total_dimers + 0.0))
        feature_row.append(math.sqrt(feature_row[1] * feature_row[2]))
        alignment = nw_hairpin.needleman_wunsch(1, 1, record[1], GtoA_sequences[record[0].split('#')[0]])
        feature_row.append(alignment[0] / (len(DNA_seq) + 0.0))
        feature_row.append(align_type[record[0]])
        feature_writer.writerow(feature_row)
        
def countKmers(DNA_seq, kmer, printout=False):
    count = 0
    for i in xrange(len(DNA_seq)):
        if DNA_seq[i:i+len(kmer)] == kmer:
            count += 1
    if printout:
        print "DNA: " + DNA_seq
        print "kmer: " + kmer
        print "count:" + count
    return count

def findEntropy(DNA_seq):
    DNA_length = len(DNA_seq)
    freq_list = [DNA_seq.count("A")/(DNA_length + 0.0), DNA_seq.count("C")/(DNA_length + 0.0), 
                  DNA_seq.count("G")/(DNA_length + 0.0), DNA_seq.count("T")/(DNA_length + 0.0),
                  DNA_seq.count("N")/(DNA_length + 0.0)]
    ent = 0
    for freq in freq_list:
        if freq != 0:
            ent += freq*math.log(freq)
    return (-ent, freq_list)
    
    
if __name__ == '__main__':
    main()