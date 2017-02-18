from Bio import SeqIO
#from Bio.SeqIO.QualityIO import PairedFastaQualIterator

fqFile = "/home/jsporter/Real_Data/SC_CHR_22_11321_CLEANER.fastq"
outfile = open("/home/jsporter/Real_Data/SC_CHR_22_11321_CLEANERER.fastq", "w")
outfile_e = open("/home/jsporter/Real_Data/SC_CHR_22_11321_BIG.fastq", "w")

cutoff = 2048

#print "Writing: "
#records = []
#num_deleted = 0
#for r in SeqIO.parse(fastq_iter, "fastq"):
#    if len(r.seq) >= 10:
#        records.append(r)
#        
#    else:
#        num_deleted += 1
#        print r
#        print ''
fastq_iter_e = SeqIO.parse(fqFile,"fastq")
records_e = (r for r in fastq_iter_e if len(r.seq) >= cutoff )# if len(r.seq) < 10)
count_e = SeqIO.write(records_e, outfile_e, "fastq")
outfile_e.close()

print "Number of records removed: " + str(count_e)

fastq_iter = SeqIO.parse(fqFile,"fastq")
records = (r for r in fastq_iter if len(r.seq) < cutoff )
count = SeqIO.write(records, outfile, "fastq")

outfile.close() 

print "Number of records merged: " + str(count)
#print "Number removed: "  + str(num_deleted)