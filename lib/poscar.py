""" Vasp poscar file.
"""
import os
from mytype import System, Molecule, Atom
from utilities import v2lattice

class Poscar():
    """ vasp poscar file
    """
    def __init__(self, filename="POSCAR"):
        self.name = ''               # name of the simulation
        self.scale = 1.0             # scale factor
        self.a = []                  # a vector
        self.b = []                  # b vector
        self.c = []                  # c vector
        self.natoms = []              # atoms numbers
        self.atomtypes = []          # atom types
        self.sd = 0                  # selective dynamics
        self.coords = []             # coordinations
        
        # Read the file
        self.read(filename)          

    def read(self, filename):
        """ read vasp poscar file
        @todo : add Direct tag
        """
        f = open(filename, 'r')
        tmp = []
        for i in f:
            if i.strip().startswith("Direct"):
                break
            elif i.strip().startswith("Selective"):
                self.sd = 1
            else:
                tmp.append(i)
        for i in f:
            self.coords.append(i.strip().split())

        self.name = tmp[0].strip()
        self.scale = float(tmp[1].strip())
        self.a = tmp[2].strip().split()
        self.b = tmp[3].strip().split()
        self.c = tmp[4].strip().split()
        if len(tmp) == 6:
            self.natoms = tmp[5].strip().split()
            self.atomtypes = self.readPotcar()
        if len(tmp) == 7:
            self.natoms = [int(j) for j in tmp[6].strip().split()]
            self.atomtypes = tmp[5].strip().split()

    def readPotcar(self, potcar="POTCAR"):
        """read atom types form POTCAR file
        """
        types = []
        f = open(potcar, 'r')
        for i in f:
            if i.strip().startswith("VRHFIN"):
                a = i.strip().split("=")[1].split(":")[0]
                types.append(a)
        f.close()
        return types

    def parser(self,):
        """ parse poscar to system
        @todo : handle coordinates acoording to Direct tag
        """
        s = System()
        scale = float(self.scale)
        a = [scale*float(i) for i in self.a]
        b = [scale*float(i) for i in self.b]
        c = [scale*float(i) for i in self.c]
        s.name = self.name
        s.pbc = v2lattice(a, b, c)
        s.natoms = self.natoms
        s.atomtypes = self.atomtypes
        s.scaleFactor = self.scale
        s.geotag = "XTLGRF 200"
        for i in range(len(self.natoms)):
            if i == 0:
                prev = 0
                now = int(self.natoms[0])
            else:
                prev = now
                now += int(self.natoms[i])
            for j in range(prev, now):
                atom = Atom()
                atom.name = self.atomtypes[i]
                x = float(self.coords[j][0])
                y = float(self.coords[j][1])
                z = float(self.coords[j][2])
                atom.xFrac = [x, y, z]
                atom.x[0] = a[0]*x + b[0]*y + c[0]*z
                atom.x[1] = a[1]*x + b[1]*y + c[1]*z
                atom.x[2] = a[2]*x + b[2]*y + c[2]*z
                if len(self.coords[j]) == 6 and self.sd == 1:
                    xr = self.coords[j][3]
                    yr = self.coords[j][4]
                    zr = self.coords[j][5]
                    if xr == "F":
                        atom.xr[0] = 1
                    if yr == "F":
                        atom.xr[1] = 1
                    if zr == "F":
                        atom.xr[2] = 1
                s.atoms.append(atom)
        return s
            
if __name__ == "__main__":
    a = Poscar()
    b = a.parser()
    
