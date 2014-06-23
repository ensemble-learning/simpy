"""
parse the rxn.log file to Rxn class, and do analysis
"""

import operator

class Rxn():
    def __init__(self, line):
        self.nstep = 0
        self.id = 0
        self.nreac = 0
        self.npro = 0
        self.reac = []
        self.pro = []
        self.reacid = []
        self.proid = []
        self.reactag = ""
        self.line = line
        self.parser()

    def parser(self,):
        tokens = self.line.strip().split(None, 3)
        self.nstep = int(tokens[0].split("_")[0])
        self.id = int(tokens[0].split("_")[1])
        self.nreac = int(tokens[1])
        self.npro = int(tokens[2])
        reac = tokens[3].split("::")[0].split("+")
        reac2 = []
        for i in reac:
            for j in i.split():
                reac2.append(j)
        pro = tokens[3].split("::")[1].split("+")
        pro2 = []
        for i in pro:
            for j in i.split():
                pro2.append(j)

        for i in range(self.nreac):
            self.reac.append(reac2[2*i])
            self.reacid.append(int(reac2[2*i + 1].strip("()")))

        for i in range(self.npro):
            self.pro.append(pro2[2*i])
            self.proid.append(int(pro2[2*i + 1].strip("()")))
        
        self.reac.sort()
        self.pro.sort()
        
        self.reactag = "_" + "_".join(self.reac)
        self.reactag += "_=_"
        self.reactag += "_".join(self.pro) + "_"
        
def parse_rxn():
    lines = []
    f = open("rxn.log", "r")
    for i in f:
        if i.startswith("#"):
            pass
        else:
            if len(i.strip()) > 0:
                lines.append(i)
    return lines

def main():
    rxns = []
    lines = parse_rxn()
    for i in lines:
        a = Rxn(i)
        rxns.append(a)

    #statics
    dist = {}
    for i in rxns:
        tag = i.reactag
        if tag in dist.keys():
            dist[tag] += 1
        else:
            dist[tag] = 1
    sorted_dist = sorted(dist.iteritems(), key=operator.itemgetter(1), reverse=True)
    fp = open("dist.log", "w")
    fp.write('\n'.join('%-80s %s' %i for i in sorted_dist))
    fp.close()
    
if __name__ == "__main__":
    main()
