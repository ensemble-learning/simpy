import os

class Reaction():
    def __init__(self, input):
        self.name = []
        self.step = 0
        self.reac = []
        self.reacid = []
        self.product = []
        self.productid = []
        self.atoms = []
        self.parser(input)
    def parser(self, input):
        tokens = input
        self.step = int(tokens[0].split("_")[0])
        self.nreac = int(tokens[1])
        self.nproduct = int(tokens[2])
        reac = tokens[3].split("::")[0].split("+")
        product= tokens[3].split("::")[1].split("+")
        for i in reac:
            item = i.strip().split()
            mol = item[0] 
            molid = int(item[1][1:-1])
            self.reac.append(mol)
            self.reacid.append(molid)
        for i in product:
            item = i.strip().split()
            mol = item[0] 
            molid = int(item[1][1:-1])
            self.product.append(mol)
            self.productid.append(molid)
    def get_atoms(self, idmap):
        vmd = 1
        for i in self.reacid:
            for j in idmap[i][2:]:
                #print self.reac[counter], idmap[i][1]
                if vmd:
                    self.atoms.append(int(j) - 1)
                else:
                    self.atoms.append(int(j))
        self.atoms.sort()

def parse_molid():
    molid = {}
    f = open("molid.out", "r")
    for i in f:
        if i.strip().startswith("#"):
            pass
        else:
            tokens = i.strip().split()
            if len(tokens) > 3:
                id = int(tokens[0])
                molid[id] = tokens[1:]
    return molid

def parse_rxn(molid):
    reactions = []
    f = open("rxn.log", "r")

    for i in f:
        if i.strip().startswith("#"):
            pass
        else:
            tokens = i.strip().split(" ",3)
            if len(tokens) > 1:
                reac = Reaction(tokens)
                reac.get_atoms(molid)
                reactions.append(reac)
    return reactions

def filter_reac(reactions, mols):

    if not os.path.exists("reactions"):
        os.mkdir("reactions")
    os.chdir("reactions")

    for i in mols:
        o = open(i + "_react.dat", "w")
        for j in reactions:
            if i in j.reac:
                line = "%8d"%j.step
                line += " | "
                line += " + ".join(j.reac)
                line += " = "
                line += " + ".join(j.product)
                line += " | "
                if j.atoms:
                    line += " ".join([str(ii) for ii in j.atoms])
                line += "\n"
                o.write(line)
        o.close()
    os.chdir("..")

def filter_product(reactions, mols):

    if not os.path.exists("reactions"):
        os.mkdir("reactions")
    os.chdir("reactions")

    for i in mols:
        o = open(i + "_product.dat", "w")
        for j in reactions:
            if i in j.product:
                line = "%8d"%j.step
                line += " | "
                line += " + ".join(j.reac)
                line += " = "
                line += " + ".join(j.product)
                line += " | "
                if j.atoms:
                    line += " ".join([str(ii) for ii in j.atoms])
                line += "\n"
                o.write(line)
        o.close()

    os.chdir("..")

def main():
    molid = parse_molid()
    reactions = parse_rxn(molid)
    mols = ["H3ON", "N2", "ON2", "H2O", "H3N", "ON", "H5O2N" ]
    filter_reac(reactions, mols)
    filter_product(reactions, mols)

if __name__ == "__main__":
    main()

