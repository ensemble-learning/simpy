def parseXvg(xvgfile):
    xvg = []
    f = open(xvgfile, 'r')
    for i in f:
        if i.strip().startswith("#"):
            pass
        elif i.strip().startswith("@"):
            pass
        else:
            if len(i.split()) >= 2:
                xvg.append(i.split())
    f.close()
    return xvg
 
def rdfIntergral(xvg):
    for i in range(len(xvg)):
        xvg[i][1] = ((float(xvg[i][1])-1)*float(xvg[i][0])*float(xvg[i][0]))
    return xvg

def outputXvg2csv(xvg, outfile):
    o = open(outfile, 'w')
    for i in xvg:
        o.write("%14s,%14.4f\n"%(i[0], i[1]))
    o.close()

if __name__ == "__main__":
    xvgfile = "nana.xvg"
    outfile = "nana.csv"
    xvg = parseXvg(xvgfile)
    xvg = rdfIntergral(xvg)
    outputXvg2csv(xvg, outfile)

    xvgfile = "nawater.xvg"
    outfile = "nawater.csv"
    xvg = parseXvg(xvgfile)
    xvg = rdfIntergral(xvg)
    outputXvg2csv(xvg, outfile)

    xvgfile = "sna.xvg"
    outfile = "sna.csv"
    xvg = parseXvg(xvgfile)
    xvg = rdfIntergral(xvg)
    outputXvg2csv(xvg, outfile)

    xvgfile = "ss.xvg"
    outfile = "ss.csv"
    xvg = parseXvg(xvgfile)
    xvg = rdfIntergral(xvg)
    outputXvg2csv(xvg, outfile)

    xvgfile = "swater.xvg"
    outfile = "swater.csv"
    xvg = parseXvg(xvgfile)
    xvg = rdfIntergral(xvg)
    outputXvg2csv(xvg, outfile)

