"""parse the integrated file into seperated single files
"""

def dumpBlock(dumpfile, outfile="dump.sep"):
    """parse the dump file into blocks
    """
    counter = 0
    tag = 0
    f = open(dumpfile, 'r')
    lines = ''
    block = []
    for i in f:
        if counter > 0 and "TIMESTEP" in i:
            o = open(outfile+"%05d"%tag+".dump", 'w')
            for j in block:
                o.write(j)
            o.close()
            block = []
            tag += 1
        block.append(i)
        counter += 1
    o = open(outfile+"%05d"%tag+".dump", 'w')
    for j in block:
        o.write(j)
    o.close()
    f.close()

def geoBlock(geofile, ext="geo"):
    """parse the geo file into blocks
    """
    f = open(geofile, 'r')
    lines = ''
    block = []
    out = ''
    counter = 0
    for i in f:
        if counter > 0:
            if "BIOGRF" in i:
                o = open("%s"%out + ".%s"%ext, 'w')
                for j in block:
                    o.write(j)
                o.close()
                block = []
            elif "XTLGRF" in i:
                o = open("%s"%out + ".%s"%ext, 'w')
                for j in block:
                    o.write(j)
                o.close()
                block = []
            elif "DESCRP" in i:
                out = i.strip().split()[-1]
            block.append(i)
        counter += 1
    o = open("%s"%out + ".%s"%ext, 'w')
    for j in block:
        o.write(j)
    o.close()

if __name__ == "__main__":
    geoBlock("geo")
