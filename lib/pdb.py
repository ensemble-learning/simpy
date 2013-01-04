""" Pdb file
"""

from mytype import System, Molecule, Atom

class Pdb():
    """ Pdb file
    """
    def __init__(self, filename):
        self.name = ''
        self.pbc = []
        self.coords = []
        self.read(filename)

    def read(self, filename):
        """ read data from pdb
        """
        f = open(filename, 'r')
        for i in f:
            if i.strip().startswith("ATOM"):
                self.coords.append(i.strip())
            elif i.strip().startswith("CRYST1"):
                self.pbc = i.strip().split()[1:7]
        f.close()
    def parser(self,):
        s = System()
        if self.name:
            s.name = self.name
        else:
            s.name = "pdbjob"
        if self.pbc:
            s.pbc = [float(i) for i in self.pbc]
        for i in self.coords:
            a = Atom()
            a.name = i[11:15].strip()
            a.x[0] = float(i[27:39])
            a.x[1] = float(i[39:47])
            a.x[2] = float(i[47:55])
            s.atoms.append(a)
        return s

if __name__ == "__main__":
    print __doc__
