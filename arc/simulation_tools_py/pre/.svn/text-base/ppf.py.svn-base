import copy
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
def outputppf(ppf, filename):
    o = open(filename, 'w')
    o.write("""#DFF:PPF
#PROTOCOL = AMBER
""")
    for i in ppf:
        if i[0].startswith('#'):
            line = ' '.join([j for j in i])
        else:
            line = ': '.join([j.strip() for j in i])+'\n'
        o.write(line)
    o.close()
    #o = open(filename[:-4] + '.eqt', 'w')
    #for (i,j) in eqt.items():
    #    line = "%s :"%i.strip() + '  '.join([m for m in eqt[i]]) +'\n'
    #    o.write(line)
    o.close()
if __name__ == "__main__":
    parameters = parsefile("RESULT.ppf")
    parm = addStar(parameters)
    outputppf(parm, "test.ppf")