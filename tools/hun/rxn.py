"""
parse the rxn.log file to Rxn class, and do analysis
@log: 
2014-06-23: fix the bug of sorting name
"""

from operator import itemgetter

class Mol():
    def __init__(self, line):
        self.id = 0
        self.natom = 0
        self.name = ''
        self.atoms = []
        self.line = line
        self.lifetime = 0
        self.parser()

    def parser(self,):
        tokens = self.line.strip().split()
        self.id = int(tokens[0])
        self.natom = int(tokens[1])
        self.name = tokens[2]
        self.atoms = [int(i) for i in tokens[3:]]

class Rxn():
    def __init__(self, line):
        self.nstep = 0
        self.id = 0
        self.nreac = 0
        self.npro = 0
        self.reac = []
        self.pro = []
        self.reacid = []
        self.proid = []
        self.reactag = ""
        self.reacidtag = ""
        self.line = line
        self.parser()

    def parser(self,):
        tokens = self.line.strip().split(None, 3)
        self.nstep = int(tokens[0].split("_")[0])
        self.id = int(tokens[0].split("_")[1])
        self.nreac = int(tokens[1])
        self.npro = int(tokens[2])
        reac = tokens[3].split("::")[0].split("+")
        reac2 = []
        for i in reac:
            for j in i.split():
                reac2.append(j)
        pro = tokens[3].split("::")[1].split("+")
        pro2 = []
        for i in pro:
            for j in i.split():
                pro2.append(j)

        for i in range(self.nreac):
            self.reac.append(reac2[2*i])
            self.reacid.append(int(reac2[2*i + 1].strip("()")))

        for i in range(self.npro):
            self.pro.append(pro2[2*i])
            self.proid.append(int(pro2[2*i + 1].strip("()")))
        
        # sor the name and id simutaniously
        self.reac, self.reacid = [list(x) for x in zip(*sorted(zip(self.reac, self.reacid), key=itemgetter(0)))]
        self.pro, self.proid = [list(x) for x in zip(*sorted(zip(self.pro, self.proid), key=itemgetter(0)))]
        
        self.reactag = "_" + "_".join(self.reac)
        self.reactag += "_=_"
        self.reactag += "_".join(self.pro) + "_"
        
        self.reacidtag = "_" + "_".join([str(i) for i in self.reacid])
        self.reacidtag += "_=_"
        self.reacidtag += "_".join([str(i) for i in self.proid]) + "_"

def parse_rxn():
    lines = []
    f = open("rxn.log", "r")
    for i in f:
        if i.startswith("#"):
            pass
        else:
            if len(i.strip()) > 0:
                lines.append(i)
    return lines

def parse_molid():
    lines = []
    f = open("molid.out", "r")
    for i in f:
        if i.startswith("#"):
            pass
        else:
            if len(i.strip()) > 0:
                lines.append(i)
    return lines
    
def output_molid_ext(mols):
    o = open("molid.ext", "w")
    for i in mols:
        o.write("%6d"%i.id)
        o.write("%6d"%i.natom)
        o.write("%20s"%i.name)
        o.write("%8d"%i.lifetime)
        o.write("\n")
    o.close()
        

def main():
    pass
    
if __name__ == "__main__":
    main()
