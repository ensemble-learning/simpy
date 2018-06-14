#!/usr/bin/env python

import os, sys

head = """PM7 +
GRAPHF AUX BONDS DENSITY PI ENPART MMOK
Jobname = mopac

"""

if len(sys.argv) > 1:
    infile = sys.argv[1]
    outfile = sys.argv[1].strip().split(".")[0] + ".mop"
    f = open(infile, "r")
    lines = f.readlines()
    f.close()
    o = open(outfile, "w")
    o.write(head)
    for i in lines[2:]:
        tokens = i.strip().split()
        if len(tokens) > 3:
            print tokens
            ele = tokens[0]
            x = float(tokens[1])
            y = float(tokens[2])
            z = float(tokens[3])
            x_flag = 1
            y_flag = 1
            z_flag = 1
            o.write("%-6s%12.4f%3d%12.4f%3d%12.4f%3d\n"
                    %(ele, x, x_flag, y, y_flag, z, z_flag))
    o.close()
    

