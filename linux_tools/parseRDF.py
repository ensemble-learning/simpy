from xvg import parseXvg

import os
import os.path

def parsefolder(folder):
    for i in os.listdir(folder):
        fullname = os.path.join(folder, i)
        if os.path.isdir(fullname):
            parsefolder(fullname)
        elif i == "sow.xvg":
            RDF[fullname] = parseXvg(fullname)

def outputRDF(rdf):
    o = open("outputrdf.csv", 'w')
    keys = rdf.keys()
    keys.sort()
    for i in keys:
        o.write(i+'\n')
        for j in rdf[i]:
            o.write("%14s,%14.4f\n"%(j[0], float(j[1])))
    o.close()

if __name__ == "__main__":
    global RDF
    RDF ={}
    parsefolder('.')
    print RDF
    outputRDF(RDF)
