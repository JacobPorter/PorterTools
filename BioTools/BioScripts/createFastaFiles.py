#!/usr/bin/python
f_in = open("/research/jsporter/Data/Karthik/1.Aligned_NotMapped_Read1AND2.txt")
CT_out = open("/research/jsporter/Data/Karthik/Aligned_NotMapped_CT.fa", 'w')
GA_out = open("/research/jsporter/Data/Karthik/Aligned_NotMapped_GA.fa", 'w')

count = 0
for line in f_in:
    line_list = line.strip("\n").split("\t")
    read_name = ">" + line_list[0]
    CT_out.write(read_name + "#" + "CT" + "\n")
    CT_out.write(line_list[1] + "\n")
    GA_out.write(read_name + "#" + "CT" + "\n")
    GA_out.write(line_list[2] + "\n")
    
print "Number of lines in file:", str(count)