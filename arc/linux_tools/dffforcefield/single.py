import os
import os.path
from exp_data_0709 import EXP_DATA

def parsePotential(filename):
    f = open(filename, 'r')
    potential = 0
    for i in f:
        if i.strip().startswith("Potential"):
            potential = i.split()[1]
    return potential

def parseMsd(filename):
    f = open(filename, 'r')
    for i in f:
        if i.strip().startswith("#"):
            pass
        elif len(i.strip().split()) == 1:
            molNo = i.strip()
            break
    return molNo

def parseFolder(folder):
    for i in os.listdir(folder):
        fullname = os.path.join(folder, i)
        if os.path.isdir(fullname):
            parseFolder(fullname)
        elif i == "potential.log":
            singleMolHv[os.path.basename(folder)] = parsePotential(fullname)
        elif i.endswith(".msd"):
            singleMolNo[os.path.basename(folder)] = parseMsd(fullname)
def parseLiquidPo(filename):
    liqduidPo = {}
    f = open(filename, 'r')
    for i in f:
        if i.strip().startswith("#"):
            pass
        else:
            a = i.strip().split(',')
            if a[-1] == "None":
                liqduidPo[a[0]] = 0
            else:
                liqduidPo[a[0]] = float(a[-1])
    return liqduidPo

def calHv(singleMolHv,liquidPo):
    o = open("hv.log", 'w')
    for i in singleMolHv.keys():
        if i in liquidPo.keys():
            print i,EXP_DATA[i][2]
            hv = (float(singleMolHv[i]) - liquidPo[i]/EXP_DATA[i][2] + 8.314*EXP_DATA[i][1]/1000)/4.184
            o.write("%30s,%10.2f,%10.2f\n"%(i, hv, EXP_DATA[i][-1]))
        else:
            print i, "not include in liquid calculation"
    o.close()

if __name__ == "__main__":
    global singleMolHv, singlMolNo
    singleMolNo ={}
    singleMolHv = {}
    liquidPo = parseLiquidPo("/home/chengtao/dffforcefield/t_amber_ver4/alkane/result.csv")
    parseFolder('.')
    calHv(singleMolHv,liquidPo)
 
