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
        self.connect = []
        self.read(filename)

    def read(self, filename):
        """ read data from pdb
        """
        if filename.endswith(".pdb"):
            self.name = filename.split(".")[0]
        else:
            self.name = filename
            
        f = open(filename, 'r')
        for i in f:
            if i.strip().startswith("ATOM"):
                self.coords.append(i.strip())
            elif i.strip().startswith("HETATM"):
                self.coords.append(i.strip())
            elif i.strip().startswith("CRYST1"):
                self.pbc = i.strip().split()[1:7]
            elif i.strip().startswith("CONECT"):
                self.connect.append(i.strip())
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
            # note need to confirm with the standary pdb format
            a.name = i[11:17].strip()
            a.x[0] = float(i[27:39])
            a.x[1] = float(i[39:47])
            a.x[2] = float(i[47:55])
            s.atoms.append(a)

        # only storage the first half of the bond matrix
        for i in self.connect:
            tokens = i.strip().split()
            a1 = int(tokens[1])
            for j in [int(k) for k in tokens[2:]]:
                a2 = j
                if a2 > a1:
                    s.connect.append([a1, a2])
        return s

if __name__ == "__main__":
    print __doc__
