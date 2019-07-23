from ase import neighborlist
from ase.io import write,read
import numpy as np
import sys

if len(sys.argv) > 1:
    atoms = read(sys.argv[1])
    atps = []

    symbols = atoms.get_chemical_symbols()
    s2n = {}
    n = 0
    for i in symbols:
        if i not in s2n.keys():
            s2n[i] = n
            n += 1
        atps.append(s2n[i])

    i = neighborlist.neighbor_list('i', atoms, 6)
    coord = np.bincount(i)

    atps_nl = []
    for i in range(len(s2n.keys())):
        atps_nl.append([])

    for i in range(len(coord)):
        ni = atps[i]
        atps_nl[ni].append(coord[i])
    for i in range(len(atps_nl)):
        print(i, max(atps_nl[i]), int(max(atps_nl[i])*1.2))
