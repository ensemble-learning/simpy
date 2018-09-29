import os
import os.path
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

def parsePPffile(ppffile):
    ppf = []
    eqt = {}
    f = open(ppffile, 'r')
    for i in f:
        if i.startswith('#'):
            ppf.append(i)
        else:
            ppf.append(i.split(':'))
    f.close()
    f = open(ppffile[:-4] + '.eqt', 'r')
    for i in f:
        if i.startswith('#'):
            pass
        elif len(i.strip()) < 1:
            pass
        else:
            eqtterm = i.split(':')[0].strip()
            eqtvalue = i.split(':')[1].split()
            eqt[eqtterm] = eqtvalue
    f.close()
    return ppf, eqt

def outputppf(ppf, eqt, filename):
    o = open(filename, 'w')
    o.write("""#DFF:PPF
#PROTOCOL = AMBER
""")
    for i in ppf:
        if i[0].startswith('#'):
            pass
        else:
            line = ': '.join([j.strip() for j in i])+'\n' 
            o.write(line)
    o.close()
    o = open(filename[:-4] + '.eqt', 'w')
    for (i,j) in eqt.items():
        line = "%s :"%i.strip() + '  '.join([m for m in eqt[i]]) +'\n'
        o.write(line)
    o.close()

def assignITP(ppf, eqt, itpfile):
    eqtM = eqt.copy()
    for i in eqtM.keys():
        eqtM[i][0] = i
    ppfM = []
    for i in ppf:
        if i[0].strip() == "N12_6":
            pass
        else:
            ppfM.append(i)
    topparm = parseItpfile(itpfile)
    for i in topparm['atomtypes']:
        term = ["N12_6", "%s"%i[0], "%8.4f,%8.4f"%((float(i[4])*11.22),(float(i[5])/4.184))]
        ppfM.append(term)
    return ppfM, eqtM

def parsefolder(folder):
    for i in os.listdir(folder):
        fullname = os.path.join(folder, i)
        if os.path.isdir(fullname):
            parsefolder(fullname)
        elif i.endswith(".itp"):
            print fullname
            ppfM, eqtM = assignITP(ppf_ref, eqt_ref, fullname)
            ppfout = os.path.join(folder, "outputAMBER.ppf")
            outputppf(ppfM, eqtM, ppfout)
if __name__ == "__main__":
    ppffile = "alcoholEther.ppf"
    global  ppf_ref, eqt_ref
    ppf_ref, eqt_ref = parsePPffile(ppffile)
    parsefolder('.')
