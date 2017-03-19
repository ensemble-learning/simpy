"""
Get surface atoms from frac_atom.dat
"""
import sys

if len(sys.argv) > 1:
    cutoff = float(sys.argv[1])
else:
    print "Using default cutoff 0.3"
    cutoff = 0.30

def read_frac_atom():
    frac = []
    f = open("frac_atom.dat", "r")
    for i in f:
        tokens = i.strip()
        if len(tokens) > 0:
            frac.append(float(tokens))
    f.close()
    return frac

def read_out_music():
    f = open("out.music", "r")
    lines = f.readlines()
    f.close()
    return lines

frac = read_frac_atom()
lines = read_out_music()


surface = []
bulk = []

for i in range(len(frac)):
    tokens = lines[i+1].strip().split()
    element = tokens[4]
    x = tokens[1]
    y = tokens[2]
    z = tokens[3]
    line = " ".join([element, x, y, z])
    line = line + '\n'
    if frac[i] < cutoff:
        bulk.append(line)
    else:
        surface.append(line)

print "find %d surface atoms and %d bulk atoms"%(len(surface), len(bulk))
o = open("surface.xyz", "w")
o.write("%d\n\n"%len(surface))
for i in surface:
    o.write(i)
o.close()

o = open("bulk.xyz", "w")
o.write("%d\n\n"%len(bulk))
for i in bulk:
    o.write(i)
o.close()
        
    
