import random, math, os, copy, time
random.seed(time.time())

f = open("POSCAR", "r")
lines = f.readlines()
f.close()
#line[0] comments
#line[1] scale factor
#line[2-4] pbc
#line[5] elements
#line[6] natoms
#line[7] Selective dynamics
#line[8] label for coordination type
#line[9] coordination

head = []
coords = []
if lines[7].startswith("Selective"): 
    for i in range(9):
        head.append(lines[i])
natoms = 0
for i in lines[6].strip().split():
    natoms += int(i)
for i in range(9,9+natoms):
    coords.append(lines[i])

x = [0.0, 1.0, 0.0, 1.0, 0.21, 0.23]
nsamples = 200
nx =  int(math.sqrt(nsamples)) + 1
dx = x[1] - x[0]
dy = x[3] - x[2]
dz = x[5] - x[4]
dx = dx/nx
dy = dy/nx

samples = random.sample(range(nx*nx), nsamples)
samples.sort()

coords_new = []
for i in samples:
    xn = i/nx
    yn = i%nx
    xnew = xn*dx + random.random()*dx
    ynew = yn*dy + random.random()*dy
    znew = x[4] + random.random()*dz 
    coords_new.append([xnew, ynew, znew])

for i in range(len(coords_new)):
    folder = "r_%03d"%i
    if not os.path.exists(folder):
        os.mkdir(folder)
    os.chdir(folder)
    o = open("POSCAR", "w")
    new_head = copy.copy(head)
    new_head[5] = new_head[5].strip() + "    H\n"
    new_head[6] = new_head[6].strip() + "    1\n"
    newline = ""
    for j in coords_new[i]:
        newline += "%18.8f"%j
    newline += " T T T\n"
    for i in new_head:
        o.write(i)
    for i in coords:
        line = i.replace("T", "F")
        o.write(line)
    o.write(newline)
    o.close()
    os.chdir("..")

