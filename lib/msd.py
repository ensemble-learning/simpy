"""transfer msd file to towhee input
@version: 1.0
@author: hawkweedcheng
@contact: chengtao@sjtu.edu.cn
"""
from mytype import Atom, Bond, Molecule, System

class Msd():
    """Basic class for msd file include (name, remarks, pdc, mol)
    functions:
    @todo: implement more functions
    """
    def __init__(self, filename):
        self.name = filename
        """@ivar: input file name
        @type: char
        """
        self.remarks = []
        """@ivar: remarks (lines starts with '#')
        @type: list
        """
        self.pbc = []
        """@ivar: pbc info
        @type: list 
        """
        self.sys = System()
        """@ivar: system
        @type: System
        """
    def __lineToAtom(self, line):
        """put each coordinate line into atom class"""
        ai = Atom()
        tokens = line.strip().split()
        ai.number = int(tokens[0])
        ai.name = tokens[1]
        ai.an = tokens[2]
        ai.type1 = tokens[3]
        ai.type2 = tokens[4]
        ai.charge = float(tokens[5])
        ai.x[0] = float(tokens[6])
        ai.x[1] = float(tokens[7])
        ai.x[2] = float(tokens[8])
        ai.resn = int(tokens[9])
        ai.resname = tokens[10]
        ai.cg = int(tokens[11])
        return ai

    def __lineToBond(self, line):
        """put each bond line into bond class
        """
        bi = Bond()
        tokens = line.strip().split()
        bi.b1 = int(tokens[0])
        bi.b2 = int(tokens[1])
        bi.type = int(tokens[2])
        return bi

    def parseFromFile(self):
        """msd file parser
        @todo: called when class is assigned
        """
        f = open(self.name, 'r')
        for i in f:
            if i.startswith('#'):
                self.remarks.append(i.strip())
            elif i.startswith('PBC'):
                for j in i.split()[1:]:
                    self.pbc.append(float(j))
            elif len(i.split()) == 12:
                self.sys.atoms.append(self.__lineToAtom(i))
            elif len(i.split()) == 3:
                self.sys.bonds.append(self.__lineToBond(i))
            else:
                pass

    def msd2Towhee(self):
        """msd file to towhee input file
        @todo: 
        """

        """unit ntype qqatom
        1    'Cl-' -1.0d0
        vibration
        0
        improper torsion
        0
        """

        o = open('fragment.towhee', 'w')
        for i in self.mol:
            o.write("unit ntype\n")
            o.write("%-5d'%5s'\n"%(i.number, i.name))
            o.write("vibration\n")
            o.write("0\n")
        o.close()

def main():
    a = Msd('sim.msd')                
    a.parseFromFile()
    a.sys.parseToMol()
    a.sys.mols[0].getNB()
    print a.sys.mols[0].nb
    

if __name__ == '__main__':
    main()
