import sys

infile = sys.argv[1] + ".gjf"
outfile = sys.argv[1] + ".xyz"

f = open(infile, "r")

n = 0

for i in f:
    if n >= 4:
        break
    n += 1
    
atoms = []
for i in f:
    tokens = i.strip().split()
    if len(tokens) == 5:
        atoms.append(tokens)
f.close()

o = open(outfile, "w")
o.write("%d\n"%(len(atoms)))
o.write("\n")

for i in atoms:
    o.write("%6s%12.4f%12.4f%12.4f\n"%(i[0],
            float(i[2]), float(i[3]), float(i[4])))
o.close()
        

