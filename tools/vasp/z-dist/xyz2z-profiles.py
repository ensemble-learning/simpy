import numpy as np

f = open("movie.xyz", "r")
lines = f.readlines()
f.close()

n_atom = int(lines[0])
n_frame = len(lines)/(n_atom + 2)

o_list = []
for i in range(n_atom+2):
    tokens = lines[i].strip().split()
    if len(tokens) == 4:
        if tokens[0] == "O":
            o_list.append(i)

o_atoms = []
for i in range(len(o_list)):
    o_atoms.append([])

for i in range(n_frame):
    for j in range(len(o_list)):
        tokens = lines[i*(n_atom + 2)+o_list[j]].strip().split()
        z = float(tokens[-1])
        o_atoms[j].append(z)

for i in range(len(o_atoms)):
    data = np.array(o_atoms[i])
    hist, bin_edges = np.histogram(data, 20)
    dx = bin_edges[1] - bin_edges[0]
    o = open("hist-o-%03d.dat"%i, "w")
    for i in range(len(hist)):
        o.write("%12.6f%8d\n"%(bin_edges[i] + 0.5*dx, hist[i]))
    o.close()
    
