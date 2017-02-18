import cPickle as Pickle

fin = open("d0L8.R1.all.ent")
ent_dict = {}

count = 0
for line in fin:
    if count <= 2:
        count += 1
        continue
    if len(line) <= 4:
        break
    spl = line.split("\t")
    zero_pos = spl[0].split("ID:")
    ent_dict[zero_pos[1].strip()] = (spl[1], spl[2], spl[3])
    count += 1
    
print count

Pickle.dump(ent_dict, open("d0L8.R1.all.ent.dict", "w"))

#   8:1101:1513:2169#0/1