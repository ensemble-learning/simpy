import os
import os.path
import copy
import shutil

def parseMSDfile(msdfile):
    f = open(msdfile, 'r')
    parm =[]
    for i in f:
        parm.append(i.strip().split())
    return parm

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
def outputMSDfile(parm, msdfile):
    shutil.copy(msdfile, msdfile+'_bak')
    o = open(msdfile, 'w')
    for i in parm:
        line = '   '.join([j for j in i])+'\n'
        o.write(line)
    o.close()

def assigncharge(msdparm, topparm):
    parm = copy.copy(msdparm)
    charges = []
    counter = 0
    for i in topparm['atoms']:
        charges.append(i[6])
    for i in parm:
        if len(i) == 12:
            i[5] = charges[counter]
            counter += 1
    return parm

def parsefolder(folder):
    for i in os.listdir(folder):
        fullname = os.path.join(folder,i)
        if os.path.isdir(fullname):
            parsefolder(fullname)
        else:
            if i.endswith(".msd"):
                msdfile = fullname
                topfile = fullname[:-4] + '.top'
                msdparm = parseMSDfile(msdfile)
                topparm = parseTopfile(topfile)
                msdparmA = assigncharge(msdparm, topparm)
                outputMSDfile(msdparmA,msdfile)
                print fullname
            else:
                pass

if __name__ == "__main__":
    parsefolder('.')
