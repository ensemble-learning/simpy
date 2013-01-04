#THis code is used to generate TATB molecules according to each molecule's center of mass

def readMolCoord(coordfile):
    a = []
    f = open( coordfile, 'r')
    for i in f:
	tokens = i.split()
        a.append([float(tokens[0]), float(tokens[1]), float(tokens[2])])
    f.close()
    return a

def calCOM(a):
    center = [0.0, 0.0, 0.0]
    for i in a:
        center[0] += i[0]
        center[1] += i[1]
        center[2] += i[2]

    center[0] /= len(a)
    center[1] /= len(a)
    center[2] /= len(a)
    return center

def generateMol(a, b, refCom):
    c = []
    for i in b:
        for j in a:
            mol = [0.0, 0.0, 0.0]
            mol[0] = j[0] - refCom[0] + i[0]    
            mol[1] = j[1] - refCom[1] + i[1]
            mol[2] = j[2] - refCom[2] + i[2]
            print mol[0]
            print mol[1]
            print mol[2]
            c.append(mol)
            print c
    return c

def writeOutfile(c):
    o = open("averageXYZ.log", 'w')
    for i in c:
        o.write('%8.3f'%i[0])
        o.write('%8.3f'%i[1])
        o.write('%8.3f'%i[2])
        o.write('\n')
    o.close()

a = readMolCoord("single.gro")
refCom = calCOM(a)
b = readMolCoord("averageCom.log")
c = generateMol(a, b, refCom)
writeOutfile(c)
