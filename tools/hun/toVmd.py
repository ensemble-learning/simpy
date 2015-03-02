f = open("molid.out", "r")
o = open("molid.vmd", "w")

"""
1 18 O6C4N8 2 106 6 41 77 42 5 109 110 105 78 114 113 1 161 165 162 166
"""

for i in f:
    if i.strip().startswith("#"):
        pass
    else:
        tokens = i.strip().split()
        n = 0
        for j in tokens:
            if n < 3:
                o.write(j+" ")
            else:
                o.write("%d "%(int(j) - 1))
            n += 1
        o.write("\n")
o.close()
f.close()
