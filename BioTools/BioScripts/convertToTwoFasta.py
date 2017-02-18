
fd = open("1.Aligned_NotMapped_Read1AND2.txt")
c_to_t = open("Karthik.ctot.fa", 'w')
g_to_a = open("Karthik.gtoa.fa", 'w')
for line in fd:
    line_list = line.strip("\n").split("\t")
    c_to_t.write(line_list[0]+"_c_to_t"+"\n")
    c_to_t.write(line_list[1]+"\n")
    g_to_a.write(line_list[0]+"_g_to_a"+"\n")
    g_to_a.write(line_list[1]+"\n")