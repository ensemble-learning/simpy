"""Gaussian 03 file
"""
import sys
from mytype import System, Molecule, Atom
from cons import ELEMENT

class G03Gjf():
    def __init__(self, filename="g03.gjf"):
        """ for gaussian 03 input file (gjf)
        """
        self.name = filename.split(".")[0]
        self.options = []
        self.methods = ''
        self.title = ''
        self.spin = 1
        self.charge = 0
        self.atoms = []
        self.connect = []
        self.redundant = []
        self.read(filename)
    def read(self, filename):
        """ read the gjf file
        """
        # if connectivity or not
        connect_flag = 0
        f = open(filename, "r")
        # read the options and methods
        for i in f:
            if len(i.strip()) > 0:
                if i.startswith("%"):
                    self.options.append(i)
                elif i.startswith("#"):
                    self.methods = i.strip()
                    if "connectivity" in i:
                        connect_flag = 1
            else:
                break

        # read the title
        for i in f:
            if len(i.strip()) > 0:
                self.title = i
            else: 
                break
        # read the spin, charge and the atomic coordinations
        for i in f:
            if len(i.strip()) > 0:
                tokens = i.strip().split()
                if len(tokens) == 2:
                    self.charge = int(tokens[0])
                    self.spin= int(tokens[1])
                else:
                    self.atoms.append(tokens)
            else:
                break
        # read connectivity if any
        if connect_flag == 1:
            for i in f:
                if len(i.strip()) > 0:
                    self.connect.append(i)
                else:
                    break
        # read redundant if any
        for i in f:
            if len(i.strip()) > 0:
                self.redundant.append(i)

        f.close()
    
    def parser(self,):
        s = System()
        s.options = self.options
        s.methods = self.methods
        for i in self.redundant:
            s.redundant.append(i.strip().split())
        s.connect = self.connect
        s.spin = self.spin
        s.charge = self.charge
        s.name = self.title.strip()
        for i in self.atoms:
            a = Atom()
            a.name = i[0]
            a.x[0] = float(i[1])
            a.x[1] = float(i[2])
            a.x[2] = float(i[3])
            s.atoms.append(a)
        return s
            
class G03LogConf():
    def __init__(self, filename='g03.log'):
        self.name = filename.split(".")[0]
        self.atoms = []
        self.stat = 0
        self.methods = []
        self.redundant = []
        self.connect = []
        self.spin = 0
        self.charge = 0
        """@var : if simulation quit normally
        """
        self.read(filename)

    def read(self, filename):
        f = open(filename, 'r')
        flag = 1
        counter = 0
        while flag:
            flag = 0
            for i in f:
                flag = 1
                if "Input orientation:" == i.strip():
                    self.atoms = []
                    counter = 0
                    break
            for i in f:
                if counter >= 4:
                    if "----" in i:
                        break
                    self.atoms.append(i.strip().split())
                counter += 1
        f.close()

        # read the keywords
        f = open(filename, 'r')
        for i in f:
            if i.strip().startswith("#"):
                self.methods = i.strip()
                break
        for i in f:
            if i.strip().startswith("-----"):
                break
        for i in f:
            if "Charge" in i and "Multiplicity" in i: 
                tokens = i.strip().split()
                self.charge = int(tokens[2])
                self.spin = int(tokens[5])
            elif "Normal termination" in i:
                self.stat = 1
            elif "The following ModRedundant input" in i:
                break

        for i in f:
            if "Iteration" in i:
                break
            else:
                self.redundant.append(i)

        for i in f:
            if "Normal termination" in i:
                self.stat = 1
        f.close()

    def __parse_keyword(self,):
        keywords = []
        tmp = ''
        state = 0
        for i in self.methods:
            if i == "#":
                pass
            else:
                if i == " ":
                    if state:
                        if len(tmp) > 0:
                            keywords.append(tmp)
                        tmp = ''
                else:
                    tmp += i
                    if i == "(":
                        state = 0
                    if i == ")":
                        state = 1
        if len(tmp) > 0:
            keywords.append(tmp)
        
        return keywords

    def parser(self,):
        s = System()
        s.name = self.name
        s.methods = self.__parse_keyword()
        s.spin = self.spin
        s.charge = self.charge
        for i in self.redundant:
            s.redundant.append(i.strip().split())
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
        sym = -1
        f = open(self.name, "r")
        counter = 0
        while(flag):
            flag = 0
            for i in f:
                flag = 1
                if r"|HF" in i:
                    sym = 0
                    lines = ""
                    lines = i.strip()
                    break
                elif r"\HF" in i:
                    sym = 1
                    lines = ""
                    lines = i.strip()
                    break
            for i in f:
                lines += i.strip()
                if len(i.strip()) == 0:
                    break
            counter += 1
        f.close()

        if sym == 0:
            for i in lines.split(r"|"):
                if "HF" in i:
                    ener = float(i[3:])
        elif sym == 1:
            for i in lines.split("\\"):
                if "HF" in i:
                    ener = float(i[3:])
        else:
            print("Error: No energy read")

        # get the zpe energy
        zpe = 0
        f = open(self.name, "r")
        for i in f:
            if "Sum of electronic and zero-point Energies" in i:
                tokens = i.strip().split("=") 
                zpe = float(tokens[1])
        f.close()
        return ener, zpe

    def getMullikenCharge(self,):
        """get Mulliken charege from QM
        """
        charges = []
        lines = []
        f = open(self.name, "r")
        for i in f:
            if "Mulliken charges" in i or \
                "Total atomic charges" in i:
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

    def getESPCharge(self,):
        """get ESP charege from QM
        """
        charges = []
        lines = []
        f = open(self.name, "r")
        for i in f:
            if "ESP charges" in i or \
                "Total atomic charges" in i:
                break
        for i in f:
            if "Sum of ESP charges" in i:
                break
            lines.append(i)
        f.close()
        for i in lines[1:]:
            tokens = i.strip().split()
            charges.append(tokens)
        return charges

def main():                
    if len(sys.argv) < 2:
        print("g03.py logfile")
    else:
        for i in sys.argv[1:]:
            geo = i.split(".")[0]
            a = G03LogConf(i)
            b = G03tools(i)
            hf, zpe = b.getEnergy()
            o = open("qm_ener.dat", "w")
            o.write("%-15s HF %.8f\n"%(geo, hf))
            o.write("%-15s ZPE %.8f\n"%(geo, zpe))
            o.close()
            charges = b.getCharge()
            o = open("qm_charge.dat", "w")
            for j in charges:
                o.write("%-15s%6s%6s%12s  # %s\n"%( geo, "1", j[0], j[2], j[1]))
            o.close()

def testff():
    a = G03Gjf() 
    a.parser()

if __name__ == "__main__":
    testff()
    main()
    
