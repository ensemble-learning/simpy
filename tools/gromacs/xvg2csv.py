#!/usr/bin/env python3

import sys

if len(sys.argv) > 1:
    fname = sys.argv[1]
    out = 'out.csv'
    d = []
    if len(sys.argv) > 2:
        out = sys.argv+'.csv'
    f = open(fname, 'r')
    for i in f:
        if i.strip().startswith('#'):
            pass
        elif i.strip().startswith('@'):
            pass
        else:
            tokens = i.strip().split()
            d.append(tokens)
    f.close()
    o = open(out, 'w')
    for i in d:
        line = ','.join(i)+'\n'
        o.write(line)
    o.close()
else: 
    print("Usage: convert xvg file to csv format")
