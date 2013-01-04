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
 
def rdfIntergral(xvg):
    for i in range(len(xvg[0])):
        xvg[0][i] = float(xvg[0][i])
        xvg[1][i] = ((float(xvg[1][i])-1)*float(xvg[0][i])*float(xvg[0][i]))
    return xvg

def aveIntergral(xvg):
    integral = []
    for i in range(len(xvg[0])):
        if xvg[0][i] > 0.85 and xvg[0][i] < 1.2:
            integral.append(scipy.integrate.simps(xvg[1][0:i],xvg[0][0:i]))
    sum = 0
    for i in integral:
        sum += i
    return sum/len(integral)
            
def outputXvg2csv(xvg, outfile):
    o = open(outfile, 'w')
    for i in range(len(xvg[0])):
        o.write("%14.4f,%14.4f\n"%(xvg[0][i], xvg[1][i]))
    o.close()

if __name__ == "__main__":
    kbunit = 7564.954981
    xvglist = ["nana", "ss", "sna", "swater", "nawater"]
    kirkwood = {}
    for i in xvglist:
        xvgfile = i + ".xvg"
        outfile = i + ".csv"
        rdfoutfile = i +"_rdf.csv"
        xvg = parseXvg(xvgfile)
        outputXvg2csv(xvg, rdfoutfile)
        xvg = rdfIntergral(xvg)
        kirkwood[i] = aveIntergral(xvg)
        outputXvg2csv(xvg, outfile)
    gcc = (kirkwood['nana'] + kirkwood['sna']*2 + kirkwood['ss'])/4*kbunit
    gcw = (kirkwood['swater'] + kirkwood['nawater'])/2*kbunit
    print "gcc    ", gcc
    print "gcw    ", gcw
    
    
