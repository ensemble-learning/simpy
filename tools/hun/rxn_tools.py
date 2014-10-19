"""
rxn tools
@todo: 
2014-06-23: add a command line parser to select different functions
@log:
2014-06-23: first version
"""

import os
import operator
#import pygraphviz as pgv
from rxn import Rxn, Mol, parse_rxn, parse_molid
from rxn import output_molid_ext

class RxnAnalysis():
    def __init__(self,):
        self.rxns = []
        self.get_rxns()

    def get_rxns(self,):
        lines = parse_rxn()
        for i in lines:
            a = Rxn(i)
            self.rxns.append(a)
    
    def get_statics(self,):
        """
        statics of the reactions
        """
        dist = {}
        for i in self.rxns:
            tag = i.reactag
            if tag in dist.keys():
                dist[tag] += 1
            else:
                dist[tag] = 1
        sorted_dist = sorted(dist.iteritems(), key=operator.itemgetter(1), reverse=True)
        fp = open("dist.log", "w")
        fp.write('\n'.join('%-80s %s' %i for i in sorted_dist))
        fp.close()

    def get_single_rxn(self, molname, nstart=-1):

        nsum = nstart
        o = open(molname + "_details.dat", "w")
        for i in self.rxns:
            n = 0 # record 
            tag = 0
            if molname in i.reac:
                n += -1 
                tag += 1
            if molname in i.pro:
                n += 1
                tag += 1
                
            if tag > 0:
                nsum += n
                o.write("%10d%50s%4d%6d\n"%(i.nstep, i.reactag, n, nsum))
        o.close()

    def get_single_rxn_filter(self, molname, nstart=-1):

        rxns = []
        nsum = nstart

        # check out the reactions involves molname
        for i in self.rxns:
            if (molname in i.reac) or (molname in i.pro):
                rxns.append(i)

        stats = [0]*len(rxns)
        pos = []
        for i in range(len(rxns)):
            if molname in rxns[i].reac:
                if molname in rxns[i].pro:
                    stats[i] = 1
                else:
                    ndx = i
                    pos.append([ndx, rxns[i].proid])
        
        neg = []
        for i in range(len(rxns)):
            if molname in rxns[i].pro and stats[i] == 0:
                ndx = i
                neg.append([ndx, rxns[i].reacid])
        
        for i in pos:
            print "-"*50
            n = 0
            if stats[i[0]] == 0:
                for j in neg:
                    if stats[j[0]] == 0:
                        for ii in i[1]:
                            for jj in j[1]:
                                print ii, jj
                                if ii == jj:
                                    #stats[i[0]] = 1
                                    #stats[j[0]] = 1
                                    n += 1
                            if n == len(i[1]):
                                print n, len(i[1]), i[0], j[0]
    """
                        if n == len(i[1]):
                            print n, len(i[1]), i[0], j[0]
                            stats[i[0]] = 1
                            stats[j[0]] = 1
    """

    """
        for i in pos:
            if stats[i[0]] == 0:
                print "%4d%27s%40s"%(stats[i[0]], rxns[i[0]].reacidtag, rxns[i[0]].reactag)

        print "-"*100
        for i in neg:
            if stats[i[0]] == 0:
                print "%4d%27s%40s"%(stats[i[0]], rxns[i[0]].reacidtag, rxns[i[0]].reactag)
    """

    def generate_dot():
        """
        generate the reaction map in .dot
        @note: need dot to get better plot
        """
        A = pgv.AGraph()
        for i in self.rxns:
            for j in range(i.nreac):
                for k in range(i.npro):
                    A.add_edge("%d_%s"%(i.reacid[j], i.reac[j]), "%d_%s"%(i.proid[k], i.pro[k]))
                    #print "%d_%s"%(i.reacid[j], i.reac[j]), "%d_%s"%(i.proid[k], i.pro[k])
        A.write('rxn_map.dot')
        B=pgv.AGraph('rxn_map.dot')
        B.layout()
        B.draw('rxn_map.png') 
        # get better svg
        os.system("dot -Tsvg rxn_map.dot -o rxn_map.svg")

    def get_lifetime():
        """
        Calculate the lifetime of species.
        """
        sim_time = self.rxns[-1].nstep
    
        mols = []
        lines = parse_molid()
        for i in lines:
            a = Mol(i)
            mols.append(a)

        for i in range(len(mols)):
            id = mols[i].id
            start = 0
            end = sim_time
            for j in range(len(self.rxns)):
                #print id, rxns[j].proid
                if id in self.rxns[j].proid:
                    start = self.rxns[j].nstep
                if id in self.rxns[j].reacid:
                    end = self.rxns[j].nstep
            mols[i].lifetime = end - start
        output_molid_ext(mols)
            
def main():
    rxns = RxnAnalysis()
    #rxns.get_statics()
    #rxns.get_single_rxn("O6C4N8", 16)
    rxns.get_single_rxn_filter("O6C4N8", 16)
    #generate_dot()
    #get_lifetime()
    
if __name__ == "__main__":
    main()
