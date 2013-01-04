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


