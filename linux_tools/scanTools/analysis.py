INPUT_TEMPLATE = """#DFF: Simulation Job
mol = dffpy.read_mol(r'%(MODEL_NAME)s')
ff, eqt = dffpy.read_ff(r'%(PPF_NAME)s')

eset = EnergySet()
efunc = ConstEnergyFunction(%(POTENTIAL)f  , "kJ/mol")
eset.append(efunc)

job = SimulationEngine()
action = GromacsPlayer(r'%(FOLDER_NAME)s/traj.trr', r'%(FOLDER_NAME)s/ener.edr')
job.append(action, 1)

action = DigestCohensiveEnergy(mol, ff, eqt)
job.append(action, 99)

job.run(mol, eset)"""


import os
import os.path
#from exp_data_0709 import EXP_DATA

def parsePotential(filename):
    result = None
    file = open(filename)
    for i in file:
        if i.startswith("Potential"):
            result = float(i.split()[1].strip())
    file.close()
    return result
    
def calcCohensiveEnergy(vars):
    input = INPUT_TEMPLATE % vars
    pin, pout = os.popen2("./dffjob")
    pin.write(input)
    pin.close()
    for i in pout:
        if i.startswith("Cohensive Energy"):
            file = open(os.path.join(vars["FOLDER_NAME"], "cohensiveE.log"), "w")
            file.write(i)
            file.close()
            pout.close()
            break
    
    
def parseFolder(folder):
    MODEL_NAME = None
    PPF_NAME = None
    POTENTIAL = None
    FOLDER_NAME=folder
    for i in os.listdir(folder):
        fullname = os.path.join(folder, i)
        if i.endswith(".msd"):
            MODEL_NAME = fullname
        elif i.endswith(".ppf"):
            PPF_NAME = fullname
        elif i == "potential.log":
            POTENTIAL = parsePotential(fullname)
        elif os.path.isdir(fullname):
            parseFolder(fullname)
    if MODEL_NAME != None and PPF_NAME != None and POTENTIAL != None:
        calcCohensiveEnergy(vars())

def readGmx(folder):
    for i in os.listdir(folder):
        fullname = os.path.join(folder, i)
        if i.endswith(".edr"):
            denLog = os.path.join(folder,"density.log")
            poLog = os.path.join(folder, "potential.log")
            os.system("echo Density | g_energy -f %s > %s"%(fullname, denLog))
            os.system("echo Potential | g_energy -f %s > %s"%(fullname, poLog))
            os.system("rm \#*")
        elif os.path.isdir(fullname):
            readGmx(fullname)

def parseDensity(file):
    f = open(file, 'r')
    density = 0
    for a in f:
        if a.startswith("Density"):
            density = a.split()[2]
    f.close()
    return density

def parseHv(file):
    f = open(file, 'r')
    for a in f:
        if a.startswith("Cohensive"):
            hv = a.split()[4]
    f.close()
    return hv

def collect(folder):
    for i in os.listdir(folder):
        fullname = os.path.join(folder, i)
        if i.strip() == "density.log":
            path = os.path.split(fullname)
            mol = '#'.join([m for m in path[0][2:].split('/')])
            density[mol] = parseDensity(fullname)
        elif i.strip() == "cohensiveE.log":
            path = os.path.split(fullname)
            mol = '#'.join([m for m in path[0][2:].split('/')])
            hv[mol] = parseHv(fullname)
        elif os.path.isdir(fullname):
            collect(fullname)
    return density, hv

def output(density, hv):
    o = open("result.csv", 'w')
    o.write("mol,density(cal),density(exp),hv(cal),hv(exp)\n")

    keys = density.keys()
    keys.sort()
    for i in keys:
        den = float(density[i])/1000
        o.write("%s,%.3f,%s\n"%(i,den,hv[i]))
        #o.write("%s,%.3f,%s,%s,%s\n"%(i,den,EXP_DATA[i][-2],hv[i], EXP_DATA[i][-1]))
    o.close()

            
if __name__ == '__main__':
    global density, hv
    density ={}
    hv = {}
    #readGmx(".")
    #parseFolder(".")
    collect(".")
    print hv.keys()
    output(density, hv)
    

