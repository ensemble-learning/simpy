"""
<<<<<<< HEAD
"""

def read_poscars():
=======
Note: require NWRITE = 1
"""

def read_poscar():
>>>>>>> 3865e352eb4a9a3803abac2f3f38a16ee714f5a8
    #get number of atoms
    f = open("POSCAR", "r")
    lines = f.readlines()
    tokens = [i for i in lines[5].strip().split()]
    atom_types = tokens
    tokens = [int(i) for i in lines[6].strip().split()]
    n_atoms = tokens
    return atom_types, n_atoms

<<<<<<< HEAD
def read_outcar():
    atom_types, atom_numbers = read_poscars()
    n_atoms = sum(atom_numbers)

    # write type.raw
    tokens = atom_numbers
    o = open("../type.raw", "w")
=======
def read_outcar(atom_numbers):
    n_atoms = sum(atom_numbers)
    box = []
    virial = []
    energy = []
    coord = []
    force = []

    f = open("OUTCAR", "r")
    flag = 1
    while(flag):
        flag = 0
        for i in f:
            if "aborting loop because EDIFF is reached" in i:
                break

        #read virial
        tokens = []
        for i in f:
            if "FORCE on cell =-STRESS in cart. coord.  units" in i:
                break

        for i in f:
            if "Total" in i:
                tokens = i.strip().split()
                break

        if len(tokens) > 0:
            xx = tokens[1]
            yy = tokens[2]
            zz = tokens[3]
            xy = yx = tokens[4]
            yz = zy = tokens[5]
            xz = zx = tokens[6]
            virial.append([xx, xy, xz, yx, yy, yz, zx, zy, zz])

        # read box
        tokens = []
        for i in f:
            if "direct lattice vectors" in i:
                break
        n = 0
        for i in f:
            n+=1
            for j in range(3):
                tokens.append(i.strip().split()[j])
            if n > 2:
                break
        if len(tokens) > 0:
            box.append(tokens)

        # read force and coords
        t1, t2 = [], []
        for i in f:
            if "TOTAL-FORCE" in i:
                break
        n = 0
        for i in f:
            if n > 0:
                t1.append(" ".join(i.strip().split()[:3]))
                t2.append(" ".join(i.strip().split()[3:]))
            n += 1
            if n >= n_atoms+1:
                break
        if len(t1) > 0:
            coord.append(t1)
            force.append(t2)

        # read energy
        tokens = []
        for i in f:
            if "FREE ENERGIE OF THE ION-ELECTRON SYSTEM" in i:
                break

        for i in f:
            flag = 1
            if "energy  without entropy" in i:
                tokens = i.strip().split()
                break
        if len(tokens) > 0:
            energy.append(tokens[-1])
    #print(len(energy), len(box), len(coord), len(force), len(virial))
    return virial, box, coord, force, energy

def write_type_rawatom(atom_numbers):
    n_atoms = sum(atom_numbers)
    tokens = atom_numbers
    o = open("type.raw", "w")
>>>>>>> 3865e352eb4a9a3803abac2f3f38a16ee714f5a8
    n = 0
    for ii in tokens:
        for jj in range(ii): 
            o.write("%d "%n)
        n += 1
    o.write("\n")

<<<<<<< HEAD
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
=======
def output_files(virial, box, coord, force, energy):
    o = open("virial.raw", "w")
    for i in virial:
        o.write(" ".join(i)+'\n')
    o.close()
    o = open("box.raw", "w")
    for i in box:
        o.write(" ".join(i)+'\n')
    o.close()
    o = open("energy.raw", "w")
    for i in energy:
        o.write(i+'\n')
    o.close()
    o = open("force.raw", "w")
    for i in force:
        o.write(" ".join(i)+'\n')
    o.close()
    o = open("coord.raw", "w")
    for i in coord:
        o.write(" ".join(i)+'\n')
    o.close()
        

def main():
    atom_types, atom_numbers = read_poscar()
    write_type_rawatom(atom_numbers)
    virial, box, coord, force, energy = read_outcar(atom_numbers)
    output_files(virial, box, coord, force, energy)
>>>>>>> 3865e352eb4a9a3803abac2f3f38a16ee714f5a8

if __name__ == "__main__":
    main()

