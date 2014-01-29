import os
import shutil
from ffield import Ffield

class Param():
    def __init__(self, filename="param", ntype=0):
        self.params = []
        self.comments = []
        self.ntype = ntype
        self.read(filename)
    def read(self, filename):
        f = open(filename, "r")
        for i in f:
            if "!" in i:
                tokens = i.strip().split("!")
                data = tokens[0].strip().split()
                comment = tokens[1]
            else:
                data = i.strip().split()
                comment = ""
            if len(data) == 6:
                self.params.append(data)
                self.comments.append(comment)
    def update(self,):
        ff = Ffield("ffield", self.ntype)
        for i in self.params:
            nsec = int(i[0]) 
            if nsec == 1:
                pass
            elif nsec == 2:
                np = int(i[1]) - 1
                nval = int(i[2])
                val = float(ff.atom[np][nval])
            elif nsec == 3:
                np = int(i[1]) - 1
                nval = int(i[2]) + 1
                val = float(ff.bond[np][nval])
            elif nsec == 4:
                np = int(i[1]) - 1
                nval = int(i[2]) + 1
                val = float(ff.off[np][nval])
            step = val*0.05
            if val > 0:
                start = val * 0.8
                end = val * 1.2
            else:
                start = val * 1.2
                end = val * 0.8
            i[3] = step
            i[4] = start
            i[5] = end
    
    def output_param(self,):
        if os.path.exists("params"):
            shutil.copy("params", "params.bak")
        o = open("params", "w")
        for i in self.params:
            o.write("%6s"%i[0])
            o.write("%6s"%i[1])
            o.write("%6s"%i[2])
            o.write("%10.4f"%i[3])
            o.write("%10.4f"%i[4])
            o.write("%10.4f\n"%i[5])
        o.close()
    
    def filter(self, a1, a2, eq):
        for i in self.comments:
            print i
          
def test():
    a = Param("params")
    a.update()
    a.output_param()

if __name__ == "__main__":
    test()
