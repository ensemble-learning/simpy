import copy
import os
import os.path
import sys
def usage():
    print "blablabla"
def parseTopfile(topfile):
    parms = {}
    term =''
    forcefield = ''
    f = open(topfile, 'r')
    for i in f:
        if i.startswith(';'):
            pass
        elif i.startswith('#'):
            forcefield = i.split()[1].strip('"')
        elif len(i.strip()) < 1:
            pass
        elif i.startswith('['):
            term = i.split()[1]
            parms[term] = []
        else:
            if term !='':
                parms[term].append(i.split())
    f.close()
    parms["forcefield"] = forcefield
    return parms

def outputTopfile(parm, filename):
    o = open(filename, 'w')
    o.write('#include "%s"\n\n'%parm["forcefield"])
    o.write('[ moleculetype ]\n')
    for i in parm['moleculetype']:
        line = ' '.join(['%10s'%a for a in i]) +'\n'
        o.write(line)
    o.write('\n[ atoms ]\n')
    for i in parm['atoms']:
        line = ' '.join(['%10s'%a for a in i]) +'\n'
        o.write(line)
    o.write('\n[ bonds ]\n')
    for i in parm['bonds']:
        line = ' '.join(['%10s'%a for a in i]) +'\n'
        o.write(line)
    o.write('\n[ pairs ]\n')
    for i in parm['pairs']:
        line = ' '.join(['%10s'%a for a in i]) +'\n'
        o.write(line)
    o.write('\n[ angles ]\n')
    for i in parm['angles']:
        line = ' '.join(['%10s'%a for a in i]) +'\n'
        o.write(line)
    o.write('\n[ dihedrals ]\n')
    for i in parm['dihedrals']:
        line = ' '.join(['%10s'%a for a in i]) +'\n'
        o.write(line)
    o.write('\n[ system ]\n')
    for i in parm['system']:
        line = ' '.join(['%10s'%a for a in i]) +'\n'
        o.write(line)
    o.write('\n[ molecules ]\n')
    for i in parm['molecules']:
        line = ' '.join(['%10s'%a for a in i]) +'\n'
        o.write(line)

    o.close()
def scaleCharge(parm,scale):
    parmM = {}
    for (i,j) in parm.items():
        parmM[i] = parm[i]
    parmM['atoms'] = []
    for j in parm['atoms']:
        m = copy.copy(j)
        m[6] = "%.4f"%(float(j[6])*scale/100)
        parmM['atoms'].append(m)
    print parmM['atoms'][1][6]
    print parm['atoms'][1][6]
    return parmM

if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()
    else:
        filename = sys.argv[1]
        runFolder = os.getcwd()
        parm = parseTopfile(filename)
        for i in [90, 100, 110]:
            if os.path.isdir("charge_%03d"%i):
                pass
            else:
                os.system("mkdir charge_%03d"%i)
            fullname = os.path.join(runFolder, "charge_%03d"%i)
            fullname = os.path.join(fullname, filename)
            parmM = scaleCharge(parm,i)
            outputTopfile(parmM, fullname)
