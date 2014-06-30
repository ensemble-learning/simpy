"""
Find the atom coordinations from the snapshots
@note: need molid.out, rxn.log and movie.xyz
@todo: move the atoms to fit the pbc.
"""

import os

class Rxn():
    def __init__(self,):
        self.tag = ''
        self.timestep = 0
        self.molid = []
        self.atomlist = []

def get_reactions(rxn):
    f = open("rxn.log", "r")
    for i in f:
        if i.strip().startswith(rxn.tag):
            tokens = i
    return tokens

def parse_reaction(rxn, tokens):
    tokens = tokens.strip().split(None, 3)
    n1 = int(tokens[1])
    n2 = int(tokens[2])
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
    for i in range(n1):
        rxn.molid.append(reac2[2*i + 1].strip("()"))

    for i in range(n2):
        rxn.molid.append(pro2[2*i + 1].strip("()"))
    
def get_atomlist(rxn):
    f = open("molid.out", "r")
    for i in f:
        tokens = i.strip().split()
        if tokens[0] in rxn.molid:
            rxn.atomlist.extend([int(j) for j in tokens[3:]])
            

def get_snapshot(rxn, folder, nframe):
    lines = []
    conf = "output%05d.xyz"%nframe
    f = open(conf, "r")
    counter = 1
    for i in f:
        if (counter - 2) in rxn.atomlist:
            lines.append(i)
        counter += 1
    f.close()
    o = open("%s/rxn_%05d.xyz"%(folder, nframe), "w")
    o.write("%d\n"%len(lines))
    o.write("\n")
    for i in lines:
        o.write(i)
    o.close()
    
def convert2pdb():
    for i in [j for j in os.listdir(".") if j.endswith("xyz")]:
        fname = i.split(".")[0]
        os.system("babel -ixyz %s.xyz -opdb %s.pdb"%(fname, fname))
    pbc = "CRYST1   12.050   10.000   16.984  90.00  92.834  90.00 P 1\n"
    for i in [j for j in os.listdir(".") if j.endswith("pdb")]:
        f = open(i, "r")
        lines = f.readlines()
        f.close()
        o = open(i, "w")
        o.write(pbc)
        for j in lines:
            o.write(j)
            o.write(j)
        o.close()
        
def main():
    rxn = Rxn()
    rxn.tag = "12788_1"
    rxn.timestep = int(rxn.tag.split("_")[0])
    tokens = get_reactions(rxn)
    parse_reaction(rxn, tokens)
    get_atomlist(rxn)
    next = 10
    folder = "rxn_%s"%rxn.tag
    if not os.path.exists(folder):
        os.mkdir("rxn_%s"%rxn.tag)
    for i in range( rxn.timestep - next,  rxn.timestep + next):
        get_snapshot(rxn, folder, i)
    os.chdir(folder)
    convert2pdb()

if __name__ == "__main__":
    main()
