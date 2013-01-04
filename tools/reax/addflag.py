""" add number flag sequentially to geo files
"""

import sys
import os
import shutil

def usage():
    print """python addflag.py flag [geo file name]
    if no geo file name, use geo as default
    """
if len(sys.argv) < 2 or len(sys.argv) > 3:
    usage()
else:
    flag = sys.argv[1]
    geofile = "geo"
    
    if len(sys.argv) == 3:
        flag = sys.argv[1]
        geofile = sys.argv[2]

    lines = []
    f = open(geofile, "r")
    counter = 0
    for i in f:
        if i.strip().startswith("DESCRP"):
            line = "DESCRP %s_%02d\n"%(flag, counter)
            counter += 1
        else:
            line = i
        lines.append(line)
        
    f.close()
    
    shutil.copy(geofile, "geo.bak")
    o = open(geofile, "w")
    for i in lines:
        o.write(i)
    o.close()
    

