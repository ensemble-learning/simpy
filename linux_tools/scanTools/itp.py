import sys
import copy
import os
import time

def usage():
    print """xxxxx
xxxxxxxxxxxxxxxxxxx
xxxxxxxxxxxxxxxxxxxxxx
"""

def parseItpfile(itpfile):
    parms = {}
    term =''
    f = open(itpfile, 'r')
    for i in f:
        if i.startswith('#') or i.startswith(';'):
            pass
        elif len(i.strip()) < 1:
            pass
        elif i.startswith('['):
            term = i.split()[1]
            parms[term] = []
        else:
            if term !='':
                parms[term].append(i.split())
    f.close()
    return parms

def writeItp(itpfile, parms):
    o = open(itpfile, 'w')
    ISOTIMEFORMAT='%Y-%m-%d %X'
    currentTime = time.strftime( ISOTIMEFORMAT, time.localtime( time.time() ) )
    o.write("#CREATED AT%s\n\n"%currentTime)
    o.write("[ defaults ]\n")
    o.write("; nbfunc   comb-rule   gen-pairs   fudgeLJ   fudgeQQ\n")
    line = "".join(["%8s"%a.strip() for a in parms['defaults'][0]])
    o.write(line+'\n')
    
    o.write("[ atomtypes ]\n")
    o.write("; name        mass      charge   ptype         sigma           epsilon\n")  
    for i in parms['atomtypes']:
        line = "".join(["%10s"%a.strip() for a in i]) + '\n'
        o.write(line)
    
    o.write("[ bondtypes ]\n")
    o.write(";   i    j func        b0            kb\n")
    for i in parms['bondtypes']:
        line = "".join(["%15s"%a.strip() for a in i]) + '\n'
        o.write(line)
    
    o.write("[ angletypes ]\n")
    o.write(";   i    j func        b0            kb\n")
    for i in parms['angletypes']:
        line = "".join(["%10s"%a.strip() for a in i]) + '\n'
        o.write(line)

    o.write("[ dihedraltypes ]\n")
    o.write(";   i    j func        b0            kb\n")
    for i in parms['dihedraltypes']:
        line = "".join(["%10s"%a.strip() for a in i]) + '\n'
        o.write(line)

    o.close()

def parmModify(parm, deltaS,epsilonS, term, atomtype):
    parmM = {}
    for (i,j) in parm.items():
        parmM[i] = parm[i]
    parmM[term] = []
    for j in parm[term]:
        m = copy.copy(j)
        if j[0] == atomtype:
            m[4] = "%.4f"%(float(j[4])*deltaS/100)
            m[5] = "%.4f"%(float(j[5])*epsilonS/100)
        parmM[term].append(m)
    return parmM
            
def scanAll(parm, itpfile):
    folder = os.getcwd()
    term = "atomtypes"
    atom = []
    for i in parm[term]:
        atom.append(i[0])
    delta = [105, 100, 95]
    epsilon = [120, 100, 80]
    for a in atom:
        for i in delta:
            for j in epsilon:
                fullname = os.path.join(folder,'%s_%03d_%03d'%(a,i,j))
                parmM = parmModify(parm, i,j, term, a)
                if os.path.isdir(fullname):
                    pass
                else:
                    os.system("mkdir %s"%fullname)
                os.system("cp dffjob.dfi *.msd run.gro run.mdp *.top grompp.sh %s"%fullname)
                outitpfile = os.path.join(fullname, itpfile)
                writeItp(outitpfile, parmM)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()
    else:
        itpfile = sys.argv[1]
        parmOr = parseItpfile(itpfile)
        scanAll(parmOr, itpfile)
