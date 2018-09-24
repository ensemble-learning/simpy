#!/usr/bin/env python

from ase.io import read, write
import sys

outfile = "output.png"

if len(sys.argv) > 1:
    infile = sys.argv[1]
    if len(sys.argv) > 2:
        outfile = sys.argv[2]
    data = read(infile)
    write(outfile, data, rotation="-60z, 120x")
