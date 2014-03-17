"""parse the ndx (gromacs) file to groups
"""

import re

class Group():
    def __init__(self, filename="index.ndx" ):
        self.name = filename
        self.names = []
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
                    self.groups[counter] += [int(j) - 1 for j in tokens]

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

def test():
    a = Group("index.ndx")
    print a.name
    print a.names
    print a.groups

if __name__ == "__main__":
    test()
