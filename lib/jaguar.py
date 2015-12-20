"""Jaguar file
"""

import sys
from mytype import System, Molecule, Atom
from cons import ELEMENT

class Jaguar():
    def __init__(self, filename="jag.in"):
        """For jaguar input file
        """
        self.name = ""
        self.options = []
        self.methods = ''
        self.title = ''
        self.spin = 1
        self.charge = 0
        self.atoms = []
        self.connect = []
        self.redundant = []
        self.read(filename)
    def read(self, filename):
        """ read .in file
        """
        f = open(filename, "r")
        for i in f:
            if i.strip().startswith("&gen"):
                break
        for i in f:
            tokens = i.strip().split()
            if i.strip().startswith("&"):
                break
        for i in f:
            if i.strip().startswith("&zmat"):
                break
        for i in f:
            tokens = i.strip().split()
            if i.strip().startswith("&"):
                break
            else:
                self.atoms.append(tokens)
            
        f.close()
    
    def parser(self,):
        s = System()
        for i in self.atoms:
            a = Atom()
            a.name = i[0]
            a.x[0] = float(i[1].strip("#"))
            a.x[1] = float(i[2].strip("#"))
            a.x[2] = float(i[3].strip("#"))
            s.atoms.append(a)
        return s
