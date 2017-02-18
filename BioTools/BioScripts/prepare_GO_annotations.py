GO_annot = open("/home/jsporter/gene_association.goa_human", "r")
GO_processed = open("/home/jsporter/gene_association.goa_human.reduced", "w")

count = 0
evid_lst = ['TAS', 'IEA', 'IDA', 'IBA', 'IMP', 'ISS', 'NAS', 'IC', 'IEP', 'IGI', 'IPI', 'EXP', 'RCA']
my_n = "\t"
GO_processed.write('orf\tgoid\thierarchy\tevidencecode\tannotation type\n')# + '\t' * 12 + '\n' )
GO_processed.flush()
for line in GO_annot:
#     if count >= 10:
#         break
    if line[0] == "!":
        pass
    else:
        my_line = line.split("\t")
#        print my_line
        if my_line[8] != "P" or my_line[6] not in evid_lst:
            continue 
        my_str = my_line[2] + my_n + my_line[4][4:] + my_n + "p" + my_n + my_line[6] + my_n + "0"
        my_str += "\t" * 12 + "\n"
        GO_processed.write(my_str)
        GO_processed.flush()
#         print my_str
        count += 1

print "Number of annotations written: %s" % str(count)

GO_annot.close()
GO_processed.close()