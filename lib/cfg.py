"""cfg file
"""

from mytype import System, Molecule, Atom

class Cfg(filename):
    """cfg file
    """
    def __init__(self, filename):
        self.name = ""
        self.pbc = []
        self.a = []
        self.b = []
        self.c = []
        self.read(filename)
        self.natoms = 0
    
    def read(self, filename):
        """ read data from cfg file
        """
        if filename.endswith(".cfg"):
            self.name = filename[:-4]
        else:
            self.name = filename
        
        f = open(filename, "r")
        self.natoms = int(f.readline().strip().split("=")[-1])
        f.readline() # unit
        for i in range(3):
            tokens = i.strip().split("=")[-1].split()[0]
            self.a.append(float(tokens))
        for i in range(3):
            tokens = i.strip().split("=")[-1].split()[0]
            self.b.append(float(tokens))
        for i in range(3):
            tokens = i.strip().split("=")[-1].split()[0]
            self.c.append(float(tokens))
        f.readline() # not sure the meaning of the line
        n = int(f.readline().strip().split("=")[-1])
        for i in range(n -3):
            f.readline()
        
        counter = 0
        for i in range(self.natoms):
            if counter % 3 == 0:
                tokens = []
            tokens.append(i.strip().split())
            self.coords.append(token)
            counter += 1

    def parser(self,):
        """ parse cfg to system
        """
        s = System()
        s.pbc = v2lattice(a, b, c)
        s.natoms = self.natoms
        for i in range(s.natoms):
            atom = Atom()
            atom.name = self.coords[i][1][0]
            x = float(self.coords[i][2][0])
            y = float(self.coords[i][2][1])
            z = float(self.coords[i][2][2])
            atom.x[0] = a[0]*x + b[0]*y + c[0]*z
            atom.x[1] = a[1]*x + b[1]*y + c[1]*z
            atom.x[2] = a[2]*x + b[2]*y + c[2]*z
            s.atoms.append(atom)

