""" DFTB gen file format
"""

SlaterKirkwood = {"O": [0.56, 3.8, 3.15], "N": [1.03, 3.8, 2.82], "C": [1.382, 3.8, 2.5], "H":[0.386, 3.5, 0.80]}

class Gen():
    def __init__(self, fname="tmp"):
        self.name = ""
        self.n = 0
        self.pbc = []
        self.atom_maps = []
        self.coords = []
        self.reader(fname,)
    def reader(self, fname):

        pbc_tag = 0
        counter = 0 
        f = open(fname, "r")
        for i in f:
            tokens = i.strip().split()
            if counter == 0:
                self.n = int(tokens[0])
                pbc_tag = tokens[1]
            elif counter == 1:
                self.atom_maps = tokens
                break
            counter += 1
        counter = 0
        for i in f:
            tokens = i.strip().split()
            self.coords.append(tokens)
            counter += 1
            if counter ==  self.n:
                break
  
        if pbc_tag == "S" or pbc_tag == "F":
            counter = 0
            for i in f:
                tokens = i.strip().split()
                if counter >0 :
                    self.pbc.append(tokens)
                if counter == 4:
                    break
                counter += 1
        f.close()
    def toSlaterKirkwood(self,):
        for i in self.coords:
            atp = int(i[1]) - 1
            params = SlaterKirkwood[self.atom_maps[atp]]
            print " ".join(["%6.2f"%j for j in params])

def test():
    import os
    os.chdir("/net/hulk/home6/chengtao/soft/simpy/testcode/gen")
    a = Gen("tkx.gen")
    a.toSlaterKirkwood()

if __name__ == "__main__":
    test()
