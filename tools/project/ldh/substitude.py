"""Substitude 1/8 Ca into Al
The starting model (model.pdb) is a 2:1 (Ca:Al) LDH structure 
from experiment with all Al converted to Ca. Then
1. replace 9 of Ca in 72 to Al (each layer).
2. exchange Al and Ca to eliminate any Al-O-Al
"""
import random

# PBC specified to model.pdb
# @improve: a variable read from pdb input

PBC = [17.324, 39.708, 33.947]

class Atom():
    def __init__(self,):
        # atom coordinations
        self.x =[0.0, 0.0, 0.0]
        # atom id
        self.id = 0

def get_dis(a1, a2, atoms):
    """ calculate the distance of two atoms
    @note: Here we can ONLY take care of orthogonal boxes.
    """
    a = PBC[0]/2
    b = PBC[1]/2
    c = PBC[2]/2
    #@note: here the atom label is a little tricky
    r1 = atoms[a1-1].x
    r2 = atoms[a2-1].x
    dx = r1[0] - r2[0]
    # take care of PBC
    dx = dx - int(dx/a)*a
    dy = r1[1] - r2[1]
    dy = dy - int(dy/b)*b
    dz = r1[2] - r2[2]
    dz = dz - int(dz/c)*c
    return dx*dx + dy*dy + dz*dz

def cal_connect(a1, sub, layer, all):
    """calculate the connectivity according to distance
    """
    cutoff = 3.6
    n = 0
    for i in range(len(sub)):
        if i == a1:
            pass
        else:
            rsq = get_dis(sub[a1], sub[i], all)
            if rsq < cutoff*cutoff:
                n += 1
    return n
    
def move(a1, sub, layer, all):
    """Exchange Ca and Al. Accept this move if result less Al-O-Al.
    Otherwise, refuse this move.
    """
    # first calculate current Al-O-Al number
    n_old = cal_connect(a1, sub, layer, all)
    x_old = sub[a1]
    # if no Al-O-Al, no necessary for any move
    if n_old > 0:
        flag = 1
        while(flag):
            flag = 0
            x_new = random.sample(layer, 1)[0]
            if x_new in sub:
                flag = 1
        sub[a1] = x_new
        n_new = cal_connect(a1, sub, layer, all)
        
        if n_new < n_old:
            return n_new
        else:
            sub[a1] = x_old
            return n_old
    else:
        return 0

def getEnergy():
    return 1

def mcmove(a1, sub, layer, all):
    """ MC move
    """
    # first calculate current Al-O-Al number
    nold = cal_connect(a1, sub, layer, all)
    assert n = 0
    x_old = sub[a1]
    ener_old = getEnergy()

    nnow = 1
    # if no Al-O-Al, no necessary for any move
    while(n):
        flag = 1
        while(flag):
            flag = 0
            x_new = random.sample(layer, 1)[0]
            if x_new in sub:
                flag = 1
        sub[a1] = x_new
        nnow = cal_connect(a1, sub, layer, all)
    ener_new = getEnergy()
    if ener_new > ener_old:
        sub[a1] = x_old

def moves(sub, layer, all, max_iter=1000):
    """make exchange moves
    @param sub: The substitued atoms (Al in this case).
    @param layer: Atoms in selected layer. 
    @param all: All of the Atoms. 
    @param max_iter: maxium iterations.
    """

    for n in range(max_iter):
        res = 0
        for i in range(len(sub)):
            tmp = move(i, sub, layer, all)
            res += tmp
        if res == 0:
            print "Step %d Res converged to 0"%n
            break
        else:
            pass
            print "Step %d Res = %d"%(n, res)

def sub(id):
    """substitude 1/8 Ca to Al and exchange Al with Ca
    to avoid any Al-O-Al bond
    """
    # In model.pdb, there are three layers L1, L2 and L3
    L1 = []
    L2 = []
    L3 = []
    atoms = []
    f = open("model.pdb", "r")
    for i in f:
        if i.startswith("ATOM"):
            atom = Atom()
            tokens = i.strip().split()
            na = int(tokens[1])
            atn = tokens[2]
            x = float(tokens[5])
            y = float(tokens[6])
            z = float(tokens[7])
            atom.id = na - 1
            atom.x[0] = x
            atom.x[1] = y
            atom.x[2] = z
            atoms.append(atom)
            if "Ca" in atn:
                if z > 7 and z < 12:
                    L1.append(na)
                elif z > 15 and z < 19:
                    L2.append(na)
                elif z > 22 and z < 27:
                    L3.append(na)
    f.close()

    #print "Ca in Layer1 %d"%len(L1)
    #print "Ca in Layer2 %d"%len(L2)
    #print "Ca in Layer3 %d"%len(L3)

    al_L1 = random.sample(L1, 9)
    al_L2 = random.sample(L2, 9)
    al_L3 = random.sample(L3, 9)
    
    moves(al_L1, L1, atoms)
    moves(al_L2, L2, atoms)
    moves(al_L3, L3, atoms)
    al = al_L1 + al_L2 + al_L3
    al.sort()
    
    f = open("model.pdb", "r")
    o = open("caldh_7_1_al_%03d.pdb"%id, "w")
    for i in f:
        line = i
        if i.strip().startswith("ATOM"):
            tokens = i.strip().split()
            na = int(tokens[1])
            if na in al:
                line = i.replace("Ca", "Al")
        o.write(line)
    o.close()

if __name__ == "__main__":
    for i in range(100):
        sub(i)
