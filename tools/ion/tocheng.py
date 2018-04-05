#!/usr/bin/env python

import os
import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", nargs="+", help="geo file name")
    parser.add_argument("-o", default="/home/tao/Dropbox", help="convert the file to other formats (geo, xyz, gjf, lammps)")
    args = parser.parse_args()
    files = args.fname
    folder = args.o
    ip = "131.215.26.219"
    usr = "tao"
    cmd = "scp -r "
    cmd += " ".join(files)
    cmd += " "+usr+"@"+ip+":"+folder
    os.system(cmd)

if __name__ == "__main__":
    main()

