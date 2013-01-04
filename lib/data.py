"""Lammps data file
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
            atom = Atom()
            atom.name = self.atomtypes[int(i[1]) - 1]
            atom.x[0] = float(i[3]) 
            atom.x[1] = float(i[4]) 
            atom.x[2] = float(i[5]) 
            s.atoms.append(atom)
        return s

if __name__ == "__main__":
    testfile = "e_2_data.data"
    a = ReaxData(testfile)