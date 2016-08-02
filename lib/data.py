"""Lammps data file
@doing: implement a general Data class
@note: remove the ReaxData after finishing Data Class
"""
import os
from mytype import System, Atom
from utilities import v2lattice
from cons import MASS2ELMENT

class ReaxData():
    """only work for reax lammps data generated from restart2data
    """
    def __init__(self, datafile="data_lammps"):
        self.name = ""
        self.natom = 0
        self.a = []
        self.b = []
        self.c = []
        self.pbc = []
        self.coords = []
        self.atomtypes = []
        self.read(datafile)
    def read(self, datafile):
        f = open(datafile, "r")
        tilt = [0, 0, 0]
        for i in f:
            tokens = i.strip().split()
            if "atoms" in i:
                # read n atoms
                self.natom = int(tokens[0])
            elif "xlo" in i:
                x = float(tokens[1]) - float(tokens[0])
            elif "ylo" in i:
                y = float(tokens[1]) - float(tokens[0])
            elif "zlo" in i:
                z = float(tokens[1]) - float(tokens[0]) 
            elif "xy" in i:        
                tilt = [float(i) for i in tokens[:3]]
            elif "Masses" in i:
                #start reading the mass
                break
        # assign a, b, c
        self.a = [x, 0, 0]
        self.b = [tilt[0], y, 0]
        self.c = [tilt[1], tilt[2], z]
        
        for i in f:
            tokens = i.strip().split()
            if len(tokens) == 0:
                pass
            elif "Atoms" in i:
                break
            else:
                self.atomtypes.append(MASS2ELMENT[int(float(tokens[1]) + 0.5)])
        
        counter = 0
        for i in f:
            if counter == self.natom:
                break
            else:
                tokens = i.strip().split()
                if len(tokens) == 0:
                    pass
                else:
                    self.coords.append(tokens)
                    counter += 1
        f.close()
        #sort the atoms in data
        tmp = []
        for i in self.coords:
            line = "%05d_"%int(i[0]) + "_".join(i)
            tmp.append(line)
        tmp.sort()
        for i in range(len(self.coords)):
            self.coords[i] = tmp[i].split("_")[1:]
        
    def parser(self):
        s = System()
        s.pbc = v2lattice(self.a, self.b, self.c)
        for i in self.coords:
            # n is used to distinguish the type of data file
            n = 0
            for j in i:
                # ignore the comments
                if "#" in j:
                    break
                n += 1
            atom = Atom()
            atom.an = int(i[0])
            if n == 7:
                # n = 7 is full type
                atom.name = self.atomtypes[int(i[2]) - 1]
                atom.x[0] = float(i[4]) 
                atom.x[1] = float(i[5]) 
                atom.x[2] = float(i[6]) 
                s.atoms.append(atom)
            elif n == 6 or n == 9:
                # n = 5 is a charge type, n = 8 is a charge type generated
                # by restart2data
                atom.name = self.atomtypes[int(i[1]) - 1]
                atom.x[0] = float(i[3]) 
                atom.x[1] = float(i[4]) 
                atom.x[2] = float(i[5]) 
                s.atoms.append(atom)

        return s

class FullData():
    """only work for reax lammps data generated from restart2data
    """
    def __init__(self, datafile="data_lammps"):
        self.name = ""
        self.n_atoms = 0
        self.n_bonds = 0
        self.n_angles = 0
        self.n_dihedrals = 0
        self.n_impropers = 0
        self.n_atomtypes = 0
        self.n_bondtypes = 0
        self.n_angletypes = 0
        self.n_dihedraltypes = 0
        self.n_impropertypes = 0
        self.a = []
        self.b = []
        self.c = []
        self.pbc = []
        self.atomtypes = []
        self.coords = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.ffparams = []
        self.elements = []
        self.read(datafile)
    def read(self, datafile):
        f = open(datafile, "r")
        tilt = [0, 0, 0]
        for i in f:
            tokens = i.strip().split()
            if "atoms" in i:
                self.n_atoms = int(tokens[0])
            elif "bonds" in i:
                self.n_bonds = int(tokens[0])
            elif "angles" in i:
                self.n_angles = int(tokens[0])
            elif "dihedrals" in i:
                self.n_dihedrals = int(tokens[0])
            elif "impropers" in i:
                self.n_impropers = int(tokens[0])
            elif "atom types" in i:
                self.n_atomtypes = int(tokens[0])
            elif "bond types" in i:
                self.n_bondtypes = int(tokens[0])
            elif "angle types" in i:
                self.n_angletypes = int(tokens[0])
            elif "dihedral types" in i:
                self.n_dihedraltypes = int(tokens[0])
            elif "improper types" in i:
                self.n_impropertypes = int(tokens[0])
            elif "xlo" in i:
                x = float(tokens[1]) - float(tokens[0])
            elif "ylo" in i:
                y = float(tokens[1]) - float(tokens[0])
            elif "zlo" in i:
                z = float(tokens[1]) - float(tokens[0]) 
            elif "xy" in i:        
                tilt = [float(i) for i in tokens[:3]]
            elif "Masses" in i:
                break

        # assign a, b, c
        self.a = [x, 0, 0]
        self.b = [tilt[0], y, 0]
        self.c = [tilt[1], tilt[2], z]
        
        for i in f:
            tokens = i.strip().split()
            if len(tokens) == 0:
                pass
            elif "Atoms" in i:
                break
            else:
                self.elements.append(MASS2ELMENT[float(tokens[1])])
                self.atomtypes.append(MASS2ELMENT[float(tokens[1])])
        
        counter = 0
        for i in f:
            if counter == self.n_atoms:
                break
            else:
                tokens = i.strip().split()
                if len(tokens) == 0:
                    pass
                else:
                    self.coords.append(tokens)
                    counter += 1

        for i in f:
            if i.strip().startswith("Bonds"):
                break

        counter = 0
        for i in f:
            if counter == self.n_bonds:
                break
            else:
                tokens = i.strip().split()
                if len(tokens) == 0:
                    pass
                else:
                    self.bonds.append([int(k) for k in tokens])
                    counter += 1

        for i in f:
            if i.strip().startswith("Angles"):
                break

        counter = 0
        for i in f:
            if counter == self.n_angles:
                break
            else:
                tokens = i.strip().split()
                if len(tokens) == 0:
                    pass
                else:
                    self.angles.append([int(k) for k in tokens])
                    counter += 1

        for i in f:
            if i.strip().startswith("Dihedrals"):
                break

        counter = 0
        for i in f:
            if counter == self.n_dihedrals:
                break
            else:
                tokens = i.strip().split()
                if len(tokens) == 0:
                    pass
                else:
                    self.dihedrals.append([int(k) for k in tokens])
                    counter += 1
        for i in f:
            if i.strip().startswith("Pair Coeffs"):
                self.ffparams.append(i)
                break
        for i in f:
            self.ffparams.append(i)

        f.close()
        """
        print self.coords, "\n\n\n"
        print self.bonds, "\n\n\n"
        print self.angles, "\n\n\n"
        print self.dihedrals, "\n\n\n"
        print self.ffparams, "\n\n\n"
        """

        #sort the atoms in data
        tmp = []
        for i in self.coords:
            line = "%05d_"%int(i[0]) + "_".join(i)
            tmp.append(line)
        tmp.sort()
        for i in range(len(self.coords)):
            self.coords[i] = tmp[i].split("_")[1:]
        
    def parser(self):
        s = System()
        s.pbc = v2lattice(self.a, self.b, self.c)
        s.atomtypes = self.atomtypes
        s.n_atoms = self.n_atoms
        s.n_bonds = self.n_bonds
        s.n_angles = self.n_angles
        s.n_dihedrals = self.n_dihedrals
        s.n_impropers = self.n_impropers
        s.n_atomtypes = self.n_atomtypes
        s.n_bondtypes = self.n_bondtypes
        s.n_angletypes = self.n_angletypes
        s.n_dihedraltypes = self.n_dihedraltypes
        s.n_impropertypes = self.n_impropertypes
        s.bonds = self.bonds
        s.angles = self.angles
        s.dihedrals = self.dihedrals
        s.ffparams = self.ffparams

        for i in self.coords:
            # n is used to distinguish the type of data file
            n = 0
            for j in i:
                # ignore the comments
                if "#" in j:
                    break
                n += 1
            atom = Atom()
            atom.an = int(i[0])
            if n == 7:
                # n = 7 is full type
                atom.resn = int(i[1])
                atom.name = self.elements[int(i[2]) - 1]
                atom.element = self.elements[int(i[2]) - 1]
                atom.type2 = int(i[2])
                atom.charge = float(i[3])
                atom.x[0] = float(i[4]) 
                atom.x[1] = float(i[5]) 
                atom.x[2] = float(i[6]) 
                s.atoms.append(atom)
        return s

def test():
    testfile = "lammps.data"
    a = FullData(testfile)

if __name__ == "__main__":
    test()
