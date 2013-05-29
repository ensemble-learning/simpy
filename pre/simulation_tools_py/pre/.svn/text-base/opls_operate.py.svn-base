def parseTop(filename):
    f = open(filename, 'r')
    pre = []
    parm = []
    end = []
    for i in f:
        if i.startswith("[ atoms ]"):
            break
        else:
            pre.append(i)
    for i in f:
        if len(i.strip()) == 0:
            break
        else:
            parm.append(i)
    for i in f:
        end.append(i)
    f.close()
    return pre, parm, end

def outputTop(pre, parm, end):
    o = open("outputOpls.top", 'w')
    for i in pre:
        o.write(i)
    o.write("[ atoms ]\n")
    for i in parm:
        o.write(i)
    o.write("\n")
    for i in end:
        o.write(i)
    o.close()

def dff2opls(parm):
    from opls import DFF2OPLS
    outParm = []
    for i in parm:
        a = i.strip().split()
        a[1] = DFF2OPLS[a[1]][0]
        a[-2] = DFF2OPLS[a[1]][1]
        outParm.append("%6d%13s%7d%13s%13s%7d%13.4f%13.5f\n"%(j for j in a))
    return outParm

if __name__ == "__main__":
    pre, parm, end= parseTop("toulene.top")
    dff2opls(parm)
    outputTop(pre, parm, end)
