#!/usr/bin/python
fd_in = open("E14-d0.L8.alignment")
fd_recover = open("E14-d0.L8.alignment.recover.fa", 'w')
fd_left = open("E14-d0.L8.alignment.left.fa", 'w')
fd_right =  open("E14-d0.L8.alignment.right.fa", 'w')

count = 0
for line in fd_in:
    if line[0] == "#":
        continue
    else:
        spl = line.split("\t")
        if float(spl[4].strip()) == 1:
            rec_id = spl[0].strip()
            fd_recover.write(">%s#r\n" % rec_id)
            fd_recover.write("%s\n" % spl[1].strip())
            fd_recover.flush()
            fd_left.write(">%s#0/1\n" % rec_id)
            fd_left.write("%s\n" % spl[2].strip())
            fd_left.flush()
            fd_right.write(">%s#0/2\n" % rec_id)
            fd_right.write("%s\n" % spl[3].strip())
            fd_right.flush()
            count += 1
            
print "Sequences written from recovery files: %d" % count 