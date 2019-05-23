""" xyz file
"""

from mytype import System, Molecule, Atom

class Xyz():
    """ xyz file
    """
    def __init__(self, filename):
        self.name = ''
        self.natoms = 0
        self.comment = ''
        self.coords = []
        self.read(filename)

    def read(self, filename):
        """ read data from pdb
        """
        f = open(filename, 'r')
        self.natoms = int(f.readline().strip())
        self.comment = f.readline()
        for i in range(self.natoms):
            self.coords.append(f.readline())
        f.close()

    def parser(self,):
        s = System()
        if self.name:
            s.name = self.name
        else:
            s.name = "xyzjob"
        for i in self.coords:
            tokens = i.strip().split()
            a = Atom()
            a.name = tokens[0]
            a.x[0] = float(tokens[1])
            a.x[1] = float(tokens[2])
            a.x[2] = float(tokens[3])
            s.atoms.append(a)
        return s

if __name__ == "__main__":
    print(__doc__)
