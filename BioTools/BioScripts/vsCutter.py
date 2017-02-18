from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
import os
import string
from os import listdir
from os.path import isfile, join

matchers = ("Strong match", "Moderate match", "Weak match", "Suspect origin")
end_file = "Length="
#no_hits = "***** No hits found *****"
start_path = "/home/jsporter/"

no_cuts = 0
begin_cuts = 0
end_cuts = 0
both_cuts = 0
discard = 0

def findContamination(lst):
    cont_list = []
    i = 0
    while string.rfind(lst[i], end_file) == -1:
        for match in matchers:
            match_ind = string.rfind(lst[i], match)
            if match_ind > -1:
                i += 1
                loc = string.split(lst[i], '\t')
                cont_list.append((int(loc[0]), int(loc[1]), match))
                j = i
                while True:
                #if match == matchers[3]:
                    j += 1
                    loc2 = string.split(lst[j], '\t')
                    try:
                        cont_list.append((int(loc2[0]), int(loc2[1]), match))
                    except Exception:
                        break
                break
        i += 1
    if (cont_list == []):
        return None
    else:
        return cont_list

def combiner(c_lst):
    result = []
    curr = list(c_lst[0])
    i = 1
    while i < len(c_lst):
        if (curr[1] + 1 == c_lst[i][0]):
            curr[1] = c_lst[i][1]
            curr[2] += " " + c_lst[i][2]
        else:
            result.append(curr)
            curr = list(c_lst[i])
        i += 1
    result.append(curr)
    return result

def cutter(c_lst, rec, fastaFile):
    global no_cuts, begin_cuts, end_cuts, both_cuts, discard
    len_rec = len(rec.seq)
    #c_lst = c_lst.sort(key=lambda tup: tup[0])
    c_lst.sort(key=lambda tup: tup[0])
    cuts = combiner(c_lst)
    if (len(cuts) >= 3):
        discard += 1
        return (None, "Discarding file: " ,cuts)
    for item in cuts:
        if item[0] != 1 and item[1] != len_rec:
            discard += 1
            return (None, "Discarding file: " ,cuts)
    if len(cuts) == 1:
        if cuts[0][0] == 1:
            begin_cuts += 1
            #print "Cutting beginning: " + fastaFile
            return (rec[cuts[0][1]:], "Cutting beginning: ", cuts)
        else:
            end_cuts += 1
            #print "Cutting ending: " + fastaFile
            return (rec[:cuts[0][0]-1], "Cutting ending: ", cuts)
    else:
        #print "Cutting both: " + fastaFile
        both_cuts += 1
        return (rec[cuts[0][1]: cuts[1][0]-1], "Cutting both: ", cuts)
    #reduce(lambda x,y: x ,c_lst, [])
    #for c1 in c_lst:


dir_list = ["001", "002", "003", "004", "005", "006", "007"]

for dr in dir_list:
    myPath = start_path + dr + "/"
    myFastaPath = myPath + "fasta/"
    myQualityPath = myPath + "qscore/"
    myVSPath = myPath + "vs/"
    myFastQPath = myPath + "fastq/"
    
    onlyfiles = [ f[:-6] for f in listdir(myFastaPath) if isfile(join(myFastaPath,f)) ]
    
    for myFl in onlyfiles:
        fastaFile = os.path.join(myFastaPath, myFl + '.fasta')
        qscoreFile = os.path.join(myQualityPath, myFl + '.qscore')
        vsFile = open(os.path.join(myVSPath, myFl + '.vs'))
        
        contam_location = findContamination(vsFile.readlines())
        records = PairedFastaQualIterator(open(fastaFile), open(qscoreFile))
        handle = open("temp.fastq", "w")
        count = SeqIO.write(records, handle, "fastq")
        handle.close()
        for rec in SeqIO.parse("temp.fastq", "fastq"):
            out = [rec, "No cuts: ", []]
            if (contam_location != None):
                out = cutter(contam_location, rec, fastaFile)
            else:
                no_cuts += 1
            break
        if out[0] != None:
            fastqFile = open(os.path.join(myFastQPath, myFl + '.fastq'), 'w')
            count = SeqIO.write(out[0], fastqFile , "fastq")
            if count != 1:
                print "Error: there can be only one sequence " + fastaFile 
        
        print out[1] + fastaFile + " " + str(out[2])

print "No cuts " +  str(no_cuts)
print "Begin cuts " + str(begin_cuts)
print "End cuts " + str(end_cuts)
print "Both cuts " + str(both_cuts)
print "Discarded " + str(discard)