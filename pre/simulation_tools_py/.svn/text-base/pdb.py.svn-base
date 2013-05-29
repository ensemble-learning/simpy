"""transfer pdb file to towhee input
@version: 1.0
@author: hawkweedcheng
@contact: chengtao@sjtu.edu.cn
@see: http://www.wwpdb.org/index.html
"""

from mol import Atom

class pdb():
    """Basic class for pdb file include (name, remarks, pdc, mol)
    functions:
    @todo: implement more functions
    """
    def __init__(self, filename):
        self.name = filename
        self.remarks = []
        self.pbc = []
        self.mol = []
    def __lineToAtom(self, line):
        """put each coordinate line into atom class"""
        ai = Atom()
        tokens = line.split()
        ai.number = int(tokens[1])
        ai.name = tokens[2]
        ai.resname = tokens[3]
        ai.x[0] = float(tokens[5])
        ai.x[1] = float(tokens[6])
        ai.x[2] = float(tokens[7])
        ai.element = tokens[-1].strip()
        ai.charge = float(tokens[-2])
        return ai

    def parseFromFile(self):
        """pdb file parser
        @todo: called when class is assigned
        """
        f = open(self.name, 'r')
        for i in f:
            if i.startswith('REMARK'):
                self.remarks.append(i.strip())
            elif i.startswith('CRYST1'):
                for j in i.split()[1:-1]:
                    self.pbc.append(float(j))
            elif i.startswith('ATOM'):
                self.mol.append(self.__lineToAtom(i))
            else:
                pass
    def pdb2Towhee(self):
        """pdb file to towhee input file
        @todo: 1. extend with connectivity 2. regulate the format
        """
        o = open('fragment.towhee', 'w')
        for i in self.mol:
            o.write("unit ntype\n")
            o.write("%-5d%7s\n"%(i.number, "'%s'"%i.name))
            o.write("vibration\n")
            o.write("0\n")
        o.close()

def main():
    a = pdb('lit.pdb')                
    a.parseFromFile()
    a.pdb2Towhee()

if __name__ == '__main__':
    main()
