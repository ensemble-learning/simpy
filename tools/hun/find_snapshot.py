"""
Find the atom coordinations from the snapshots
"""

import os, sys
import subprocess
import socket
import shutil

LIB = ''

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simpy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simupy/lib"
elif socket.gethostname() == "atom.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "ion.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant12":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"

sys.path.insert(0 , LIB)

from block import xyzBlock

class Rxn():
    def __init__(self,):
        self.tag = ''
        self.timestep = 0
        self.molid_reac = []
        self.molid_pro = []
        self.atomlist_reac = []
        self.atomlist_pro = []

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
        rxn.molid_reac.append(reac2[2*i + 1].strip("()"))

    for i in range(n2):
        rxn.molid_pro.append(pro2[2*i + 1].strip("()"))
    
def get_atomlist(rxn):
    f = open("molid.out", "r")
    for i in f:
        tokens = i.strip().split()
        if tokens[0] in rxn.molid_reac:
            rxn.atomlist_reac.extend([int(j) for j in tokens[3:]])
        if tokens[0] in rxn.molid_pro:
            rxn.atomlist_pro.extend([int(j) for j in tokens[3:]])
    # Sort the list
    rxn.atomlist_reac.sort()
    rxn.atomlist_pro.sort()

def get_snapshot(rxn, folder, nframe):
    lines = []
    conf = os.path.join(os.getcwd(), "conf")
    conf = os.path.join(conf, "output%05d.xyz"%nframe)
    f = open(conf, "r")
    counter = 1
    for i in f:
        if (counter - 2) in rxn.atomlist_reac:
            lines.append(i)
        counter += 1
    f.close()
    o = open("%s/rxn_%05d.xyz"%(folder, nframe), "w")
    o.write("%d\n"%len(lines))
    o.write("\n")
    for i in lines:
        o.write(i)
    o.close()
    
def convert2pdb(t0, dt):
    for i in range(t0 - dt, t0 + dt):
        fname = "rxn_%05d"%i
        os.popen("babel -ixyz %s.xyz -opdb %s.pdb"%(fname, fname))
    pbc = "CRYST1   11.027   16.042   10.430  71.16  88.34  87.51 P 1\n"
    for i in range(t0 - dt, t0 + dt):
        fname = "rxn_%05d"%i
        f = open(fname +".pdb", "r")
        lines = f.readlines()
        f.close()
        o = open(fname +".pdb", "w")
        o.write(pbc)
        for j in lines:
            o.write(j)
        o.close()
        #subprocess.Popen("genconf -f %s.pdb -o %s.pdb -nbox 2 2 2"%(fname, fname), \
        #                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        os.system("genconf -f %s.pdb -o %s.pdb -nbox 2 2 2"%(fname, fname))
        f.close()
        
def find_rxn(rxn_tag):
    rxn = Rxn()
    rxn.tag = rxn_tag
    rxn.timestep = int(rxn.tag.split("_")[0])
    tokens = get_reactions(rxn)
    parse_reaction(rxn, tokens)
    get_atomlist(rxn)
    next = 10
    folder = "rxn_%s"%rxn.tag
    if not os.path.exists(folder):
        os.mkdir("rxn_%s"%rxn.tag)
    for i in range( rxn.timestep - next,  rxn.timestep + next):
        #print "---------------------------working %d--------------------------"%i
        get_snapshot(rxn, folder, i)
    os.chdir(folder)
    convert2pdb(rxn.timestep, next)
    #os.system("cat rxn*.pdb > total.pdb")
    for i in rxn.atomlist_reac:
        print i-1,
    os.chdir("..")

def main():
    print "Reading the xyz file ..."
    shutil.copy("movie.xyz", "conf")
    if not os.path.exists("conf"):
        os.mkdir("conf")
    os.chdir("conf")
    n = 192
    xyzBlock("movie.xyz", n, )
    os.remove("movie.xyz")
    os.chdir("..")

    rlist = []
    f = open("rxn.log", "r")
    for i in f:
        if i.strip().startswith("#"):
            pass
        else:
            tokens = i.strip().split()
            if len(tokens) > 1:
                rxn = tokens[0]
                rlist.append(rxn)
    f.close()
    for i in rlist:
        find_rxn(i)
    
if __name__ == "__main__":
    main()
