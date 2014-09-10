"""parse the ndx (gromacs) file to groups
@log: 
Thu Jun 19 20:55:09 PDT 2014
output the file to index.ndx
"""

import re
import os
import shutil

class Group():
    def __init__(self, filename="index.ndx" ):
        self.name = filename
        self.names = []
        self.ngroups = 0
        self.groups = []
        self.subgroups = []
        self.grpnatoms = [] # number of atoms in each group
        self.read(filename)
        
    def read(self, filename):
        name = ''
        pattern = re.compile(r'\[\s*(\w+)\s*\]')
        f = open(filename, "r")
        counter = -1
        for i in f:
            if i.strip().startswith("["):
                match = pattern.match(i)
                if match:
                    name = match.group(1)
                    self.names.append(name)
                    self.groups.append([])
                    counter += 1
            else:
                tokens = i.strip().split()
                if len(tokens) > 0:
                    self.groups[counter].extend([int(j) for j in tokens])
        self.ngroups = counter + 1
        for i in range(self.ngroups):
            self.grpnatoms.append(self.groups[i])

    def tosubgroups(self,): 
        """Catlog the atoms into subgroups
        """
        for i in range(len(self.grpnatoms)):
            self.subgroups.append([])
            natom = self.grpnatoms[i]
            counter = 0
            subgrp = []
            for j in self.groups[i]:
                if counter % natom == 0:
                    if counter > 0:
                        self.subgroups[i].append(subgrp)
                        subgrp = []
                subgrp.append(j)
                counter += 1
            self.subgroups[i].append(subgrp)

def output(ndx):
    """
    output to index file
    """
    if os.path.exists("index.ndx"):
        shutil.copy("index.ndx", "index.ndx.bak")
    o = open("index2.ndx", "w")
    for i in range(ndx.ngroups):
        o.write("[ %s ]\n"%ndx.names[i])
        counter = 0
        for j in ndx.groups[i]:
            if counter > 0 and counter % 15 == 0:
                o.write("\n")
            o.write("%5d"%j)
            counter += 1
        o.write("\n")
            
    o.close()

def test():
    a = Group("index.ndx")
    print a.name
    print len(a.names), a.ngroups
    output(a)

if __name__ == "__main__":
    test()
