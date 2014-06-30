"""
seperate the packed cif files into single file.
Usage: python sep_cif.py fname
"""

import sys
import os

def main():
    """
    Seperate the cif files. Here we use "_database_code" tag to seperate the file
    """
    ciffile = sys.argv[1]

    flist = []
    data = []
    f = open(ciffile, "r")
    
    lines = []
    for i in f:
        if i.strip().startswith("data_"):
            fname = i.strip()
            flist.append(fname)
        lines.append(i)
        if i.strip().startswith("_database_code"):
            data.append(lines)
            lines = []
    f.close()
        
    for i in range(len(data)):
        o = open(flist[i] + ".cif", "w")
        for j in data[i]:
            o.write(j)
        o.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print __doc__
    else:
        main()

