""" Generate the n kpoints according to the reciprocal box
"""
import numpy as np
import sys
from ase.io import read

ACC = 30

def get_kpoints(fname="POSCAR"):
    atoms = read(fname)
    cell = atoms.get_reciprocal_cell()

    k = []
    for i in range(3):
        a = np.linalg.norm(cell[i])
        print a
        a = int(np.ceil(a * ACC))
        if a%2 == 1:
            a = a + 1
        k.append(a)

    o = open("KPOINTS", "w")
    o.write("""K-Points
 0
Monkhorst Pack
%d  %d  %d
 0  0  0
"""%(k[0], k[1], k[2]))
    o.close()

def main():
    if len(sys.argv) > 1:
        get_kpoints(sys.argv[1])
    else:
        get_kpoints()

if __name__ == "__main__":
    main()
    
    



