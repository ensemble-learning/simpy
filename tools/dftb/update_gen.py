"""update the coordinations from input to dftb_in.hsd
"""

import sys
import shutil

def get_hsd():
    f = open("dftb_in.hsd", "r")

    head = []
    end = []

    for i in f:
        head.append(i)
        if i.strip().startswith("Geometry"):
            break

    for i in f:
        if i.strip().startswith("}"):
            end.append(i)
            break

    for i in f:
        end.append(i)

    f.close()
    return head, end

def get_update(fname):
    geo = []
    f = open(fname, "r")
    for i in f:
        if len(i.strip()) > 0:
            geo.append(i)
    f.close()
    return geo

def main():
    fname = sys.argv[1]
    geo = get_update(fname)
    head, end = get_hsd()
    lines = head + geo + end
    
    shutil.copy("dftb_in.hsd", "dftb_in.hsd.bak")
    o = open("dftb_in.hsd", "w")
    for i in lines:
        o.write(i)
    o.close()

if __name__ == "__main__":
    main()
