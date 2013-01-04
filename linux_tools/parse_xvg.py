#!/usr/bin/env python

def parse_xvg(fname):
    x = []
    y = []
    f = open(fname, 'r')

    for i in f:
        if len(i.strip()) == 0:
            pass
        elif i.strip().startswith('@'):
            pass
        elif i.strip().startswith('#'):
            pass
        else:
            tokens = i.split()
            x.append(float(tokens[0]))
            y.append(float(tokens[-1]))
    
    f.close()
    return x, y

if __name__ == "__main__":
    print "test"
