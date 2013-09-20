""" Sort the atoms according the connetivity in .mdf
1. Get the atom index. The atoms are distinguish to 
Ca and C (CO3), and save into two lists. For CO3, the
sequence of atoms are [c, o, o, o, c, o, o, o....]
2. Read the original sequence from .car.
3. Output the sorted coordinations
@ref: 09072013
@ver: 1.01
@log: 09182013: extend the mdf parser
@log: 09192013: to simpy and renamed as car.py
"""
from mytype import System, Molecule, Atom

class Atom_mdf():
    def __init__(self,):
        self.name = ""
        self.element = ""
        self.connect = []
        
class Mdf():
    def __init__(self, mdffile):
        self.atoms = []
        self.atom_names = []
        self.read(mdffile)
    def read(self, mdffile):
        atoms = []
        f = open(mdffile, "r")
        for i in f:
            if i.strip().startswith("@"):
                pass
            elif i.strip().startswith("!"):
                pass
            elif i.strip().startswith("#"):
                pass
            else:
                tokens = i.strip().split()
                # 11 is kind of abitrary here.
                if len(tokens) > 11:
                    atoms.append(tokens)
        self.parser(atoms)
    def parser(self, atoms):
        for i in atoms:
            atom = Atom_mdf()
            atom.name = i[0]
            atom.element = i[1]
            self.atom_names.append(atom.name)
            if len(i) > 12:
                for j in i[12:]:
                    atom.connect.append(j)
            self.atoms.append(atom)
                
class Car():
    def __init__(self, carfile):
        self.pbc_tag = 0
        self.pbc = []
        self.coords = []
        self.read( carfile)
    def read(self, carfile):

        f = open(carfile, "r")
        counter = 0
        #!BIOSYM archive 3
        #PBC=ON
        #Materials Studio Generated CAR File
        #!DATE Tue Sep 17 23:14:47 2013
        #PBC   26.1570   26.1570   25.2780   90.0000   90.0000   90.0000 (P1)
        for i in f:
            if counter == 0:
                pass
            elif counter == 1:
                if "PBC=ON" in i:
                    self.tag_pbc = 1
            elif counter == 2:
                pass
            elif counter == 3:
                pass
            elif counter == 4:
                self.pbc = i
            else:
                tokens = i.strip().split()
                # here we do not consider "end" tag
                # 5 is not accurate
                if len(tokens) > 5:
                    self.coords.append(i)
            counter += 1
    def parser(self,):
        s = System()
        if self.pbc:
            s.pbc = [float(i) for i in self.pbc.split()[1:7]]
        for i in self.coords:
            a = Atom()
            a.name = i[:6].strip()
            a.x[0] = float(i[6:20])
            a.x[1] = float(i[20:35])
            a.x[2] = float(i[35:50])
            s.atoms.append(a)
        return s
            
def test():
    pass
    
if __name__ == "__main__":
    test()
