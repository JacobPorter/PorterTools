#!/usr/bin/python
import optparse, sys, os
import csv

#from os import listdir
#from os.path import isfile, join


def main():
    p = optparse.OptionParser()
    p.add_option('--input','-i',help='Directory of bootstrap sample alignment results files.')
    options,_ = p.parse_args()
#    if len(args) == 0:
#        print "Use --help for help documentation"
#        return
    print "Computing bootstrap p-value with files in directory: %s" % options.input
    sys.stdout.flush()
    computePValue(options.input, 0.05)
    print "Done!"
    
def getSampleNumber(fasta_dir, file_name):
    return file_name.split(".")[2]
    

def scanAlignmentOutput(fd):
    stat_dict = {}
    first_boolean = True
    for line in fd:
        if line.startswith("Sequences analysed in total:"):
            number = int(line.split(":")[1].strip())
            if first_boolean:
                stat_dict["CtoT_total"] = number
            else:
                stat_dict["GtoA_total"] = number
        elif line.startswith("Number of alignments with a unique best hit from the different alignments:"):
            number = int(line.split(":")[1].strip())
            if first_boolean:
                stat_dict["CtoT_unique"] = number
                first_boolean = False
            else:
                stat_dict["GtoA_unique"] = number
        elif line.startswith("Total reads aligned:"):
            number = int(line.split(":")[1].strip())
            stat_dict["Recovered_total"] = number
        elif line.startswith("Unique reads:"):
            number = int(line.split(":")[1].strip())
            stat_dict["Recovered_unique"] = number
    return stat_dict

def computePValue(fasta_dir, cutoff):
    onlyfiles = [ f for f in os.listdir(fasta_dir) if os.path.isfile(os.path.join(fasta_dir,f)) ]
    
    header = ["Sample_Num", "Total_Sequences_Analysed", "CtoT_Bismark_Alignment", "GtoA_Bismark_Alignment", "Average_Bisulfite_Alignment", "Recovered_Bowtie2_Alignment", "Is_Recovered_Worse?"]
    f_str = os.path.join(fasta_dir, "bootstrapAlignmentResults.csv")
    tsvfile = open(f_str, 'w')
    tsvwriter = csv.writer(tsvfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
    tsvwriter.writerow(header)
    sig_count = 0
    sample_count = 0
    for f in onlyfiles:
        sample_num = getSampleNumber(fasta_dir, f)
        stat_dict = scanAlignmentOutput(open(os.path.join(fasta_dir, f)))
        if not (stat_dict["CtoT_total"] == stat_dict["GtoA_total"] and stat_dict["GtoA_total"] == stat_dict["Recovered_total"]):
            print "Check the calculations.  File %s does not have equal total alignment counts."
        total_seqs = stat_dict["CtoT_total"]
        CtoT_align_perc = stat_dict["CtoT_unique"] / (total_seqs + 0.0)
        GtoA_align_perc = stat_dict["GtoA_unique"] / (total_seqs + 0.0)
        recovered_align_perc = stat_dict["Recovered_unique"] / (total_seqs + 0.0)
        avg_BS = (CtoT_align_perc + GtoA_align_perc) / (2.0)
        row_temp = [sample_num, stat_dict["CtoT_total"], CtoT_align_perc, GtoA_align_perc, avg_BS, recovered_align_perc]
        row_temp.append(0)
        this_is_sig = 0 if row_temp[4] < row_temp[5] else 1
        sig_count += this_is_sig
        row_temp[6] = this_is_sig
        tsvwriter.writerow(row_temp)
        sample_count += 1
    p_value = sig_count / (sample_count + 0.0)
    print "The total number of samples is %s and the total number of results with better bisulfite treated (unique) alignment efficiency is %s" % (str(sample_count), str(sig_count))
    print "The p-value is: %s" % str(p_value)
    if p_value >= cutoff:
        print "NOT significant!"
    else:
        print "YES significant!"
        

if __name__ == '__main__':
    main()
