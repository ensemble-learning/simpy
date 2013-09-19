import sys
import shutil

fname = sys.argv[1]

f = open(fname, "r")
o = open("tmp.car", "w")

for i in f:
    line = i
    tokens = i.strip().split()
    if len(tokens) == 9:
        atom = tokens[7]
        line = i.replace("xx", "%-2s"%atom)
    o.write(line)
        
o.close()
f.close()

shutil.copy(fname, "backup.car")
shutil.copy("tmp.car", fname)
