"""parse the integrated file into seperated single files
"""

def dumpBlock(dumpfile, outfile="dump.sep", dt = 1):
    """parse the dump file into blocks
    @return: the number of frames in dump file.
    """
    counter = 0
    tag = 0
    f = open(dumpfile, 'r')
    lines = ''
    block = []
    for i in f:
        if counter > 0 and "TIMESTEP" in i:
            if tag % dt == 0:
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
    return tag

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

def g03Block(g03file, ext="log"):
    """parse the g03 log file into blocks
    """
    head = []
    block_prev = []
    block_now = []
    f = open(g03file, 'r')

    for i in f:
        if i.strip().startswith("Berny optimization"):
            break
        head.append(i)

    block_now.append(i)

    counter = 0
    flag = 1
    while flag:
        flag = 0
        for i in f:
            flag = 1
            if i.strip().startswith("Berny optimization"):
                block_prev = block_now
                block_now = []
                break
            if i.strip().startswith("Step number   1"):
                ofile = "%s_%02d.%s"%(g03file.split(".")[0], counter, ext)
                o = open(ofile, "w")
                for j in block_prev:
                    o.write(j)
                o.close()
                counter += 1
            block_now.append(i)

    ofile = "%s_%02d.%s"%(g03file.split(".")[0], counter, ext)
    o = open(ofile, "w")
    for j in block_prev:
        o.write(j)
    o.close()

    print counter


if __name__ == "__main__":
    #g03Block("/home/tao/h3_scan_4.log")
    geoBlock("geo")

