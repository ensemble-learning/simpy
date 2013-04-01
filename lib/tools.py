""" Utilities to manipulate, sort, add or remove atoms
"""

from mytype import System, Molecule, Atom
from data import ReaxData
import random

class LdhProject():
    """ Build the crystal structure starting from CaO crystal to 
    match the experiment XRD data.
    """
    def __init__(self, system):
        self.name = system.name
        # here we only consider orthogonal box
        self.pbc = system.pbc
        self.atoms = system.atoms
        self.r = 2.4076 # the distance between Ca-O in CaO crystal
        self.re_al = 8  # the number of Ca replaced into Al 
    def re1(self):
        """The first replace scheme: replace Ca to Al and H
        to obtain Ca: Al = 7: 1 LDO
        @todo: exclude Al-O-Al bond
        """
        n_al = self.re_al
        n_re = 2 * n_al
        ca_layer = []
        counter = 0
        for i in self.atoms:
            name = i.name
            natom = i.an
            x = i.x[0]
            y = i.x[1]
            z = i.x[2]
            if name == "Ca":
                nl = int((z+0.01)/self.r)
                if nl + 1 > len(ca_layer):
                    ca_layer.append([])
                ca_layer[nl].append(natom)
            counter += 1
        # define the Al and H here
        # random select replaced Ca
        # here re is the total replaced Ca (Al + H)
        re_layer = []
        al_layer = []
        for i in range(len(ca_layer)):
            al_layer.append([])
            re_layer.append([])
            
        for i in range(len(ca_layer)):
            re_layer[i] = random.sample(ca_layer[i], n_re)
            al_layer[i] = random.sample(re_layer[i], n_al)
        
        # relalce the label: First H
        for i in self.atoms:
            natom = i.an
            for j in re_layer:
                for k in j:
                    if natom == k:
                        i.name = "H"
        # replace the label (continue): Al
        for i in self.atoms:
            natom = i.an
            for j in al_layer:
                for k in j:
                    if natom == k:
                        i.name = "Al"
        # print summary
        print "The number of Ca in total: %d"%counter
        print "There are %d layers. In each layer, there are %d Ca"%(len(ca_layer), len(ca_layer[0]))
        print "%d Ca replaced into %d Al and %d H"%(n_re, n_al, n_re - n_al)
        
def test():
    """test code"""
    testfile = "a666.data"
    a = ReaxData(testfile)
    b = a.parser()
    c = LdhProject(b)
    c.re1()
    
if __name__ == "__main__":
    test()