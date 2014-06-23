"""
rxn tools
@todo: 
2014-06-23: add a command line parser to select different functions
@log:
2014-06-23: first version
"""

import operator
import pygraphviz as pgv
from rxn import Rxn, parse_rxn 

def get_statics():
    """
    statics of the reactions
    """
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

def generate_dot():
    """
    generate the reaction map in .dot
    @note: need dot to get better plot
    """
    A = pgv.AGraph()
    rxns = []
    lines = parse_rxn()
    for i in lines:
        a = Rxn(i)
        rxns.append(a)
    
    for i in rxns:
        for j in range(i.nreac):
            for k in range(i.npro):
                A.add_edge("%d_%s"%(i.reacid[j], i.reac[j]), "%d_%s"%(i.proid[k], i.pro[k]))
                #print "%d_%s"%(i.reacid[j], i.reac[j]), "%d_%s"%(i.proid[k], i.pro[k])
    A.write('rxn_map.dot')
    B=pgv.AGraph('rxn_map.dot')
    B.layout()
    B.draw('rxn_map.png') 
    # get better svg
    #os.system("dot -Tsvg rxn_map.dot -o rxn_map.svg")

def main():
    #get_statics()
    generate_dot()
    
if __name__ == "__main__":
    main()
