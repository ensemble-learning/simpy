#!/usr/bin/env python

"""build seperated sub-folders for each simulation
"""
import os
import shutil

def usage():
    print """bsr: build seperated runs
    usage:  python bsr.py
    req: need a flist file include all the datafile
    """

def main():
    f = open("flist", "r")
    for i in f:
        fname = i.strip()
        folder = fname.split(".")[0]
        if not os.path.exists(folder):
            os.mkdir(folder)
        shutil.copy(fname, folder)

if __name__ == "__main__":
    main()
