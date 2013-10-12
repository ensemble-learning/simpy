class MyForceField():
    def __init__(self,):
        self.pair_coeffs = []
        self.bond_coeffs = []
        self.angle_coeffs = []
        self.charges = []
        self.read()
    def read(self, ffield="ffield"):
        f = open(ffield, "r")
        for i in f:
            if i.strip().startswith("Pair"):
                break
        for i in f:
            if i.strip().startswith("Bond"):
                break
            self.pair_coeffs.append(self.__parse_line(i))
        for i in f:
            if i.strip().startswith("Angle"):
                break
            self.bond_coeffs.append(self.__parse_line(i))
        for i in f:
            if i.strip().startswith("Charge"):
                break
            self.angle_coeffs.append(self.__parse_line(i))
        for i in f:
            self.charges.append(self.__parse_line(i))
    def __parse_line(self, line):
        val = line.strip().split("#")[0]
        return val + "\n"
        
class LammpsData():
    def __init__(self,):
        self.sect1 = []
        self.masses = []
        self.pair_coeffs = []
        self.bond_coeffs = []
        self.angle_coeffs = []
        self.coords = []
        self.sect2 = []
        self.read()

    def read(self,):
        f = open("lammps.data", "r")
    
        self.sec1 = []
        for i in f:
            if "Masses" in i:
                break
            self.sec1.append(i)
        self.masses = []
        for i in f:
            if "Pair Coeffs" in i:
                break
            tokens = i.strip().split()
            if len(tokens) > 1:
                self.masses.append(i)
        
        self.pair_coeffs = []
        for i in f:
            if "Bond Coeffs" in i:
                break
            tokens = i.strip().split()
            if len(tokens) > 1:
                self.pair_coeffs.append(i)
        
        self.bond_coeffs = []
        for i in f:
            if "Angle Coeffs" in i:
                break
            tokens = i.strip().split()
            if len(tokens) > 1:
                self.bond_coeffs.append(i)
        
        self.angle_coeffs = []
        for i in f:
            if "Atoms" in i:
                break
            tokens = i.strip().split()
            if len(tokens) > 1:
               self.angle_coeffs.append(i)
        
        self.coords = []
        for i in f:
            if "Bonds" in i:
                break
            tokens = i.strip().split()
            if len(tokens) > 1:
                self.coords.append(i)
        
        self.sec2 = []
        for i in f:
            self.sec2.append(i)
        
        f.close()

    def assign_ff(self,):
        npair = len(self.pair_coeffs)
        nbond = len(self.bond_coeffs)
        nangle = len(self.angle_coeffs)
        natom = len(self.coords)
        ff = MyForceField()

    def assign_charge(self,):
        ff = MyForceField()
        charges = [float(i) for i in ff.charges]
        for i in range(len(self.coords)):
            tokens = self.coords[i].strip().split()
            atp =  int(tokens[2]) - 1
            tokens[0] = "%7s"%tokens[0]
            tokens[1] = "%7s"%tokens[1]
            tokens[2] = "%4s"%tokens[2]
            tokens[3] = "%11.6f"%charges[atp]
            tokens[4] = "%16.9s"%tokens[4]
            tokens[5] = "%16.9s"%tokens[5]
            tokens[6] = "%16.9s"%tokens[6]
            self.coords[i] = "".join(tokens) + "\n"
         
    def output(self,):
        o = open("output.data", "w")
        for i in self.sec1:
            o.write(i)
        o.write("\nMasses\n\n")
        for i in self.masses:
            o.write(i)
        o.write("\nPair Coeffs\n\n")
        for i in self.pair_coeffs:
            o.write(i)
        o.write("\nBond Coeffs\n\n")
        for i in self.bond_coeffs:
            o.write(i)
        o.write("\nAngle Coeffs\n\n")
        for i in self.angle_coeffs:
            o.write(i)
        o.write("\nAtoms Coeffs\n\n")
        for i in self.coords:
            o.write(i)
        o.write("\nBonds\n\n")
        for i in self.sec2:
            o.write(i)
        o.close()
    
    
if __name__ == "__main__":
    a = LammpsData()
    a.read()
    a.assign_charge()
    a.output()
        
