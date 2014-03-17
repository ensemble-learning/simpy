""" Generate the n kpoints according to the reciprocal box
"""
import numpy as np
from ase.io import read

ACC = 30

def get_kpoints():
    atoms = read("POSCAR")
    cell = atoms.get_reciprocal_cell()

    k = []
    for i in range(3):
        a = np.linalg.norm(cell[0])
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
    get_kpoints()

if __name__ == "__main__":
    main()
    
    



