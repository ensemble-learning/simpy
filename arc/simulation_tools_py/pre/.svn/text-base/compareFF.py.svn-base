# 2009-11-24
# By chengtao
# compare two set of force field parameters and check out the ones with differences
#parsefile handle the inputfile
#compareParm compare the force field parameter by force field terms
#checkParm check if the parameters are the same
# 2009-12-04
#ffversion check out force field parameters according force field version

import copy
ModelFile = "alcoholEther.ppf"
TargetFile = "AMBER20.ppf"
def parsefile(file):
    #split every line by ':'
    #ignore the lines begin with '#' or 'AT'
    parameters = []
    f = open(file)
    for i in f:
        if i.strip().startswith("#") or i.strip().startswith("AT"):
            pass
        else:
            parameters.append(i.split(":"))
    return parameters

def addStar(parameters):
    parm = copy.copy(parameters)
    for i in parm:
        list = []
        for j in i[2].split(','):
            if '*' in j.strip():
                list.append(j.strip())
            else:
                list.append(j.strip()+'*')
        i[2] = ",".join(["%8s"%m for m in list])
    return parm
                
def compareParm(t1, t2):
    if t1[0].strip() == t2[0].strip():
    #with the same function style
        counter = 0
        a1 = t1[1].split(",")
        a2 = t2[1].split(",")
        num = len(a1)
        for i in range(num):
            if a1[i].strip() == a2[i].strip():
                counter += 1
        if counter == num:
            #print t1
            #print t2
            return True
        else:
            counter = 0
        for i in range(num):
            if a1[i].strip() == a2[num-1-i].strip():
                counter += 1
        if counter == num:
            #print t1
            #print t2
            return True
        else:
            return False
    else:
        return False

def checkParm(t1, t2):
    for i in t2:
        if compareParm(t1,i) == True:
            return True
    return False

def ffversion(ff, ver):
    ffver = []
    for i in ff:
        if len(i)>=4:
            if float(i[3])>=ver:
            #check out force field parameters higher than ver or equal to ver
                ffver.append(i)
            else:
                pass
                #print i
        else:
            print "error in "
    return ffver
if __name__ == "__main__":
    M = parsefile(ModelFile)
    T = parsefile(TargetFile)
    ver = 2.0
    Tver = ffversion(T, ver)
    lines = []
    for i in M:
        if checkParm(i, Tver) == False:
            lines.append(":".join( [term for term in i]).strip() + ":3.0\n")
            #this lines repesents new parameters
        else:
            lines.append(":".join( [term for term in i]))
            #lines.append(":".join( [term for term in i]).strip() + ":2.0\n")
    o = open("compose.ppf", 'w')
    for i in lines:
        o.write(i)
    o.close()
        