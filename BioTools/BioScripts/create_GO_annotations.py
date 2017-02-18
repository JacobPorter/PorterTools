GO_annot = open("/home/jsporter/liver_data_Gene_GO_associations.tsv", "r")
GO_processed = open("/home/jsporter/liver_3vs2cell_Gene_GO_associations.processed.tsv", "w")

count = 0
#next = ""
evid_lst = ['TAS', 'IEA', 'IDA', 'IBA', 'IMP', 'ISS', 'NAS', 'IC', 'IEP', 'IGI', 'IPI', 'EXP', 'RCA']
for line in GO_annot:
    if count == 0:
        count +=1
        continue
#     if count >= 10:
#         break
    my_l = line.split("\t")
    
    if (my_l[4][0] == 'b' or my_l[4][0] == 'B') and my_l[5] in evid_lst:
        my_str = my_l[1] + "\t" + my_l[6][4:] + "\t" + my_l[7][:-2] + "\t" + "p" + "\t" + my_l[5] + "\t" + "1"
        my_str += "\n"
        #print [my_l[7][:-2]]
        #print my_str
        GO_processed.write(my_str)
        GO_processed.flush()
        count += 1
    else:
        continue
    
print "Number of associations created: " + str(count)
GO_annot.close()
GO_processed.close()