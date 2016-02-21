import numpy as np

def read_neighbour_list():
    nlist = []
    f = open("nlist.dat", "r")
    for i in f:
        tokens = i.strip().split()
        if len(tokens) > 0:
            atoms = [int(j) for j in tokens]
            nlist.append(atoms)
    return nlist
            
def cal_general_coordination(sur, nlist):
    gcn = []
    for i in range(len(sur)):
        if sur[i] == 1:
            c = []
            c_sum = 0.0
            for j in nlist[i]:
                c.append(cn[j-1])
                c_sum += cn[j-1]
            c_max = max(c)
            g = c_sum/c_max
            #print g, c_sum, c_max
            gcn.append(g)
        else:
            gcn.append(-1)
    o = open("gcn.dat", "w")
    o1 = open("gcn_sur.dat", "w")
    for i in range(len(gcn)):
        o.write("%8.2f\n"%gcn[i])
        if gcn[i] > 0:
            o1.write("%8.2f\n"%gcn[i])
    o.close()
    o1.close()
    return gcn
    
sur = np.loadtxt('sur_sas.dat')
cn = np.loadtxt('cn.dat') 
nlist = read_neighbour_list()
gcn = cal_general_coordination(sur, nlist)



