"""
"""

def read_poscars():
    #get number of atoms
    f = open("POSCAR", "r")
    lines = f.readlines()
    tokens = [i for i in lines[5].strip().split()]
    atom_types = tokens
    tokens = [int(i) for i in lines[6].strip().split()]
    n_atoms = tokens
    return atom_types, n_atoms

def read_outcar():
    atom_types, atom_numbers = read_poscars()
    n_atoms = sum(atom_numbers)

    # write type.raw
    tokens = atom_numbers
    o = open("../type.raw", "w")
    n = 0
    for ii in tokens:
        for jj in range(ii): 
            o.write("%d "%n)
        n += 1
    o.write("\n")

    f = open("OUTCAR", "r")

    # write box.raw
    tokens = []
    for i in f:
        if "direct lattice vectors" in i:
            break
    n = 0
    for i in f:
        n+=1
        tokens.append(i.strip().split())
        if n > 2:
            break

    o = open("box.raw", "w")
    for ii in tokens:
        for jj in range(3):
            o.write(ii[jj]+" ")
    o.write("\n")
    o.close()

    # write coord.raw
    tokens = []
    for i in f:
        if "position of ions in cartesian coordinates" in i:
            break

    n = 0
    for i in f:
        tokens.append(i.strip().split())
        n += 1
        if n >= n_atoms:
            break
    
    o = open("coord.raw", "w")
    for ii in tokens:
        for jj in range(3):
            o.write(ii[jj]+" ")
    o.write("\n")
    o.close()

    #write virial.raw
    tokens = []
    for i in f:
        if "FORCE on cell =-STRESS in cart. coord.  units" in i:
            break

    for i in f:
        if "Total" in i:
            tokens = i.strip().split()
            break

    o = open("virial.raw", "w")
    xx = tokens[1]
    yy = tokens[2]
    zz = tokens[3]
    xy = yx = tokens[4]
    yz = zy = tokens[5]
    xz = zx = tokens[6]
    virals = [xx, xy, xz, yx, yy, yz, zx, zy, zz]
    for ii in virals:
        o.write(ii + " ")
    o.write("\n")
    o.close()

    # write force.raw
    tokens = []
    for i in f:
        if "TOTAL-FORCE" in i:
            break
    n = 0
    for i in f:
        tokens.append(i.strip().split())
        n += 1
        if n >= n_atoms+1:
            break

    o = open("force.raw", "w")
    for ii in tokens[1:]:
        for jj in range(3,6):
            o.write(ii[jj] + " ")
    o.write("\n")
    o.close()

    # write energy.raw
    tokens = []
    for i in f:
        if "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" in i:
            break

    for i in f:
        if "energy  without entropy" in i:
            tokens = i.strip().split()
            break

    o = open("energy.raw", "w")
    o.write(tokens[-1])
    o.write("\n")
    o.close()

    f.close()

def main():
    read_outcar()

if __name__ == "__main__":
    main()

