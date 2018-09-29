import scipy.integrate

def parseXvg(xvgfile):
    xvg = [[],[]]
    f = open(xvgfile, 'r')
    for i in f:
        if i.strip().startswith("#"):
            pass
        elif i.strip().startswith("@"):
            pass
        else:
            if len(i.split()) >= 2:
                xvg[0].append(float(i.split()[0]))
                xvg[1].append(float(i.split()[1]))
    f.close()
    return xvg
 
def aveIntergral(xvg):
    integral = [[],[]]
    for i in range(len(xvg[0])):
        if xvg[0][i] > 0.85 and xvg[0][i] < 1.2:
            integral[0].append(xvg[0][i])
            integral[1].append(xvg[1][i])
    integralR = scipy.integrate.simps(integral[1],integral[0]))
    return integralR
            
def outputXvg2csv(xvg, outfile):
    o = open(outfile, 'w')
    for i in range(len(xvg[0])):
        o.write("%14.4f,%14.4f\n"%(xvg[0][i], xvg[1][i]))
    o.close()

if __name__ == "__main__":
    xvglist = ["a%02d_d"%i for i in range(1, 11)]
    trans = {}
    for i in xvglist:
        xvgfile = i + ".xvg"
        outfile = i + ".csv"
        xvg = parseXvg(xvgfile)
        outputXvg2csv(xvg, outfile)
        trans[i] = aveIntergral(xvg)
        outputXvg2csv(xvg, outfile)
    print trans
    
    
    
