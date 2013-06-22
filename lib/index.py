"""parse the ndx (gromacs) file to groups
"""

import re

class Group():
    def __init__(self, filename="index.ndx" ):
        self.name = filename
        self.names = []
        self.groups = []
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

def test():
    a = Group("index.ndx")


if __name__ == "__main__":
    test()
