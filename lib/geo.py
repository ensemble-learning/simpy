""" reaxFF geo file format parser
log:
11-10-2012: include bgf head (the first line in bgf;
"""

from mytype import System, Molecule, Atom

class Geo():
    """reaxFF geo file format
    """
    def __init__(self, filename):
        self.name = ''
        self.type = ''
        self.pbc = []
        self.coords = []
        self.read(filename)
    
    def read(self, filename):
        """ read geo file
        """
        f = open(filename, 'r')
        for i in f:
            if i.strip().startswith("XTLGRF"):
                self.type = i.strip()
            elif i.strip().startswith("BIOGRF"):
                self.type = i.strip()
            elif i.strip().startswith("DESCRP"):
                self.name = i.strip().split()[-1]
            elif i.strip().startswith("CRYSTX"):
                self.pbc = i.strip().split()[1:]
            elif i.strip().startswith("HETATM"):
                self.coords.append(i.strip().split())
            else:
                pass
    def parser(self,):
        """ parse geo file into System
        """
        s = System()
        s.name = self.name
        s.pbc = [float(i) for i in self.pbc]
        s.geotag = self.type
        for i in self.coords:
            a = Atom()
            a.name= i[2]
            a.x[0] = float(i[3])
            a.x[1] = float(i[4])
            a.x[2] = float(i[5])
            s.atoms.append(a)
        return s

def debug():
    """ debug the code
    """
    testfile = "../../debug/geo"
    a = Geo(testfile)
    b = a.parser()
    b.assignAtomTypes()

if __name__ == "__main__":
    debug()

    
