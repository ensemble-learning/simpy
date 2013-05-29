"""charge related

"""

import os
name = os.path.basename(os.getcwd())

charges = []
f = open("bonds.reax", "r")

for i in f:
    if i.strip().startswith("#"):
        pass
    else:
        if len(i.strip()) > 0:
            tokens = i.strip().split()
            charges.append(tokens[-1])
f.close()

atps = []

f = open("geo", "r")
for i in f:
    if i.strip().startswith("HE"):
        tokens = i.strip().split()
        atps.append(tokens[2])
f.close()

for i in range(len(charges)):
    if i == 0:
        print "%-12s%8s%10s"%(name, atps[i], charges[i])
    else:
        print "%-12s%8s%10s"%('', atps[i], charges[i])
        
