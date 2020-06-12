import sys
import copy
import os
import time

def usage():
    print("""xxxxx
xxxxxxxxxxxxxxxxxxx
xxxxxxxxxxxxxxxxxxxxxx
""")

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
            if not term in parms:
                parms[term] = []
        else:
            if term !='':
                clear_comments = 1
                if clear_comments:
                    if ';' in i:
                        parms[term].append(i.split(';')[0].split())
                    else:
                        parms[term].append(i.split())
    f.close()
    return parms

def writeItp(itpfile, parms):
    o = open(itpfile, 'w')
    ISOTIMEFORMAT='%Y-%m-%d %X'
    currentTime = time.strftime( ISOTIMEFORMAT, time.localtime( time.time() ) )
    #o.write("#CREATED AT%s\n\n"%currentTime)
    
    sections = ['moleculetype', 'atoms', 'bonds', 'angles', 'dihedrals']
    for i in sections:
        if i in parms.keys():
            o.write("[ %s ]\n"%i)
            for ii in parms[i]:
                line = ''
                for iii in ii:
                    line += iii + ' '
                line += '\n'
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
                os.system("cp run.gro run.mdp *.top grompp.sh %s"%fullname)
                outitpfile = os.path.join(fullname, itpfile)
                writeItp(outitpfile, parmM)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()
    else:
        itpfile = sys.argv[1]
        parmOr = parseItpfile(itpfile)
        writeItp('test.itp', parmOr)
        #scanAll(parmOr, itpfile)
