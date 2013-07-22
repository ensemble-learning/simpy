"""
Check the consistency between geo and trainset
"""

class Reac():
    """reactions in trainset ENERGY section
    """
    def __init__(self,):
        self.mols = []
        self.ratio = []
        self.nline = 0

def check_geo():
    """check the geo to make sure no dumplicated cases in geo
    """
    mols = []
    f = open("geo", "r")
    counter = 0
    for i in f:
        if i.strip().startswith("DESCRP"):
            name = i[7:].strip()
            if name in mols:
                print "line %8d %s has been included previously!"%(counter, name)
                f.close()
                mols = []
                break
            else:
                mols.append(name)
        counter += 1
    f.close()
    return mols

def parse_trainset():
    """parse the trainset ENERGY sections into reactions
    """
    f = open("trainset.in", "r")
    for i in f:
        if i.strip().startswith("ENERGY"):
            break

    reactions = []
    counter = 1
    for i in f:
        if i.strip().startswith("#"):
            pass
        elif i.strip().startswith("ENDENERGY"):
            break
        else:
            # parse the reactions here
            reac = Reac()
            IN = 0
            OUT = 0
            mol = ''
            i = i.rsplit(None, 1)[0]
            for j in i:
                if j == "+" or j == "-" and IN == 0:
                    IN = 1
                    OUT = 0
                elif j == "/":
                    OUT = 1
                    IN = 0
                    reac.mols.append(mol.strip())
                    mol = ''
                else:
                    if IN == 1 and OUT == 0:
                        mol += j
            reac.nline = counter
            reactions.append(reac)
        counter += 1
    return reactions

def check_trainset():
    """check the consistency between geo and trainset
    """
    mols = check_geo()
    geo_tags = [0]*len(mols)
    reactions = parse_trainset()
    for i in reactions:
        n = 0
        for j in i.mols:
            if j not in mols:
                n += 1
            else:
                geo_tags[mols.index(j)] = 1
        if n > 0:
            print i.nline, i.mols
    """
    for i in range(len(geo_tags)):
        if geo_tags[i] == 0:
            print mols[i]
    """

if __name__ == "__main__":
    check_trainset()
