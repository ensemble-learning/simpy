import os
import shutil
from ffield import Ffield
import argparse

class Param():
    def __init__(self, filename="param", ntype=0):
        self.params = []
        self.ffval = []
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
    def checkout(self,):
        flist = []
        for i in os.listdir("."):
            if i.startswith("params.ffield"):
                flist.append(i)
        n = 0
        if len(flist) > 0:
            flist.sort()
            n = int(flist[-1].split(".")[-1]) + 1
        o = open("params.ffield.%02d"%n, "w")
        for i in range(len(self.ffval)):
            o.write("%4s"%self.ffval[i][0])
            o.write("%4s"%self.ffval[i][1])
            o.write("%4s"%self.ffval[i][2])
            o.write("%9.4f"%self.ffval[i][3])
            o.write(" ! %s\n"%self.comments[i])
        o.close()
        
    def update(self, ff="ffield"):
        ff = Ffield(ff, self.ntype)
        scale = 0.2
        for i in self.params:
            val = 0.0
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
            item = [i[0], i[1], i[2], val]
            self.ffval.append(item)
            step = val*0.05
            if val > 0:
                start = val * (1-scale)
                end = val * (1+scale)
            else:
                start = val * (1-scale)
                end = val * (1+scale)
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
          
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="params", nargs="?", help="param file name")
    parser.add_argument("-ff", default="ffield", nargs=1, help="ffield name")
    parser.add_argument("-checkout",action="store_true" , help="check out the force field parameters")
    parser.add_argument("-update",action="store_true" , help="update the force field parameters")

    args = parser.parse_args()
    paramfile = args.fname
    
    ff = args.ff

    a = Param(paramfile, 0)
    if args.checkout:
        a.update(ff)
        a.checkout()
    if args.update:
        a.update(ff)
        a.output_param()

if __name__ == "__main__":
    main()
