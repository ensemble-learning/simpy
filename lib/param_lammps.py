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
            if len(i.strip()) > 0:
                data = []
                if "!" in i:
                    tokens = i.strip().split("!")
                    comment = tokens[1]
                    data.append(tokens[0].split(':')[0].strip())
                    data.extend(tokens[0].split(':')[1].strip().split())
                elif "#" in i:
                    tokens = i.strip().split("#")
                    comment = tokens[1]
                    data.append(tokens[0].split(':')[0].strip())
                    data.extend(tokens[0].split(':')[1].strip().split())
                else:
                    tokens = i.strip().split(':')
                    data.append(tokens[0].strip())
                    data.extend(tokens[1].strip().split())
                    comment = ""
                if len(data) >= 6:
                    self.params.append(data)
                    self.comments.append(comment)
    def checkout(self,):
        pass
        
    def update(self, ff="ffield"):
        pass
    
    def output_param(self,):
        pass
    
    def filter(self, a1, a2, eq):
        pass
          
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="params", help="param file name")
    #parser.add_argument("-ff", default="ffield", help="ffield name")
    #parser.add_argument("-checkout",action="store_true" , help="check out the force field parameters")
    #parser.add_argument("-update",action="store_true" , help="update the force field parameters")
    #parser.add_argument("-type", nargs=1, type=int, help="Force field type: 0 for vdw; 1 for lg_inner wall")

    args = parser.parse_args()
    paramfile = args.fname
    
    ff = args.ff

    if args.type:
        ff_type = args.type[0]
    else:
        ff_type = 0
    a = Param(paramfile, ff_type)
    if args.checkout:
        a.update(ff)
        a.checkout()
    if args.update:
        a.update(ff)
        a.output_param()

if __name__ == "__main__":
    main()
