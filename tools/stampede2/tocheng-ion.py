#!/usr/bin/env python

import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", nargs="+", help="geo file name")
    parser.add_argument("-o", default="/net/hulk/PMD/tcheng/download", help="convert the file to other formats (geo, xyz, gjf, lammps)")
    args = parser.parse_args()
    files = args.fname
    folder = args.o
    ip = "ion.wag.caltech.edu"
    usr = "chengtao"
    cmd = "scp -r "
    cmd += " ".join(files)
    cmd += " "+usr+"@"+ip+":"+folder
    os.system(cmd)

if __name__ == "__main__":
    main()

