"""Gaussian 03 file
"""
import sys
from mytype import System, Molecule, Atom
from cons import ELEMENT

class G03LogConf():
    def __init__(self, filename='g03.log'):
        self.atoms = []
        self.read(filename)

    def read(self, filename):
        f = open(filename, 'r')
        flag = 1
        counter = 0
        while flag:
            flag = 0
            for i in f:
                flag = 1
                if "Standard orientation:" == i.strip():
                    self.atoms = []
                    counter = 0
                    break
            for i in f:
                if counter >= 4:
                    if "-------------" in i:
                        break
                    self.atoms.append(i.strip().split())
                counter += 1
        f.close()
    def parser(self,):
        s = System()
        for i in self.atoms:
            a = Atom()
            a.name = ELEMENT[int(i[1])]
            a.x[0] = float(i[3])
            a.x[1] = float(i[4])
            a.x[2] = float(i[5])
            s.atoms.append(a)
        return s

class G03tools():
    def __init__(self, filename="g03.log"):
        self.name = filename
    def getEnergy(self, ):
        ener = 0
        f = open(self.name, "r")
        lines = ""
        for i in f:
            if r"\HF" in i:
                lines = i.strip()
        for i in f:
            lines += i.strip()
        f.close()
        
        for i in lines.split("\\"):
            if "HF" in i:
                ener = float(i[3:])
        return ener

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "g03.py logfile"
    else:
        for i in sys.argv[1:]:
            a = G03LogConf(i)
            b = G03tools(i)
            print b.getEnergy()


