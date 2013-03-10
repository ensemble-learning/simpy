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
        """get the final energy(hf) and zpe energy from g03 output file
        """
        ener = 0
        flag = 1
        f = open(self.name, "r")
        counter = 0
        while(flag):
            flag = 0
            for i in f:
                flag = 1
                if r"|HF" in i:
                    lines = ""
                    lines = i.strip()
                    break
            for i in f:
                lines += i.strip()
                if len(i.strip()) == 0:
                    break
            counter += 1
        f.close()

        for i in lines.split(r"|"):
            if "HF" in i:
                ener = float(i[3:])
        # get the zpe energy
        zpe = 0
        f = open(self.name, "r")
        for i in f:
            if "Sum of electronic and zero-point Energies" in i:
                tokens = i.strip().split("=") 
                zpe = float(tokens[1])
        f.close()
        return ener, zpe

    def getCharge(self,):
        """get Mulliken charege from QM
        """
        charges = []
        lines = []
        f = open(self.name, "r")
        for i in f:
            if "Mulliken atomic charges" in i:
                break
        for i in f:
            if "Sum of Mulliken charges" in i:
                break
            lines.append(i)
        f.close()
        for i in lines[1:]:
            tokens = i.strip().split()
            charges.append(tokens)
        return charges

                
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print "g03.py logfile"
    else:
        for i in sys.argv[1:]:
            geo = i.split(".")[0]
            a = G03LogConf(i)
            b = G03tools(i)
            hf, zpe = b.getEnergy()
            o = open("qm_ener.dat", "w")
            o.write("%-15s HF %.4f\n"%(geo, hf))
            o.write("%-15s ZPE %.4f\n"%(geo, zpe))
            o.close()
            charges = b.getCharge()
            o = open("qm_charge.dat", "w")
            for j in charges:
                o.write("%-15s%6s%6s%12s  # %s\n"%( geo, "1", j[0], j[2], j[1]))
            o.close()

            



