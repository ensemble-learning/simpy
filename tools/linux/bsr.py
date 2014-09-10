#!/usr/bin/env python

"""build seperated sub-folders for each simulation
"""
import os
import shutil
import argparse

def usage():
    print """bsr: build seperated runs
    usage:  python bsr.py
    req: need a flist file include all the datafile
    """

def main(args):
    f = open("flist", "r")
    n = 0
    for i in f:
        fname = i.strip()
        if args.r:
            folder = fname.split(".")[-1]
        elif args.n:
            folder = "case%02d"%n
        else:
            folder = fname.rsplit(".", 1)[0]
        if folder == fname:
            folder = folder + "ext"
        if not os.path.exists(folder):
            os.mkdir(folder)
        shutil.copy(fname, folder)
        n += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", action="store_true", help="revert: using ext other than name")
    parser.add_argument("-n", action="store_true", help="number: using number other than name")
    args = parser.parse_args()
    main(args)
