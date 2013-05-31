"""Generate tables of Comparison of atomic charges between ReaxFF and QM
"""

def getCharge():
    """get charges from fort.99
    """
    f = open("fort.99", "r")
    cs = []
    counter = -1 
    for i in f:
        if "Charge" in i:
            q1 = float(i[63:72])
            q2 = float(i[72:84])
            if i.strip().startswith("Charge"):
                atom = int(i[13:17]) 
            else:
                mol = i[:20].strip()
                atom = 1
                cs.append([])
                counter += 1
                cs[counter].append(mol)
            cs[counter].append([atom, q1, q2])
    f.close()
    return cs

def toLatex(cs):
    """to latex
    """
    for i in cs:
        for j in range(1, len(i)):
            print i[0], '&', i[j][0], '&', i[j][1], '&', \
            i[j][2], r'\\'

cs = getCharge()
toLatex(cs)

