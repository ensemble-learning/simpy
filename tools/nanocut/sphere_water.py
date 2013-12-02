"""
Cut the water sphere from a water box

"""

import sys

def usage():
    print "python sphere_water.py [radius]"

def main():
    if len(sys.argv) < 2:
        usage()
    else:
         rcut = float(sys.argv[1])
         select(rcut)

def select(rcut):

    cutoff = rcut * rcut
    coords = []
    head = []
    tail = []

    center = [25.0, 25.0, 25.0]

    f = open("mixture.pdb", "r")

    for i in f:
        if i.strip().startswith("ATOM"):
            coords.append(i)
            break
        else:
            head.append(i)

    for i in f:
        if i.strip().startswith("END"):
            tail.append(i)
            break
        else:
            coords.append(i)

    f.close()

    selected = []

    n = 0
    for i in range(0, len(coords), 3):
        tokens = coords[i].strip().split()
        x = float(tokens[5]) 
        y = float(tokens[6]) 
        z = float(tokens[7]) 
        dx = x - center[0]
        dy = y - center[0]
        dz = z - center[0]
        r2 = dx*dx + dy*dy + dz*dz
        if r2 < cutoff:
            selected.append(i)
            selected.append(i + 1)
            selected.append(i + 2)
            n += 1

    print n

    o = open("sphere_%02d.pdb"%(int(rcut)), "w")

    for i in head:
        o.write(i)

    for i in selected:
        o.write(coords[i])

    for i in tail:
        o.write(i)

    o.close()

if __name__ == "__main__":
    main()
