#!/usr/bin/env python

import sys
import math
import getopt
from parse_xvg import parse_xvg

def usage():
    """Calculate the error of data according to 
    block average.
    The format is :
    block.py -t filetype -i inputfile -o outputfile

    -h : help 
    -t : input filetype;
    -i : input file
    -o : output file
    
    Ref.
    H. Flyvbjerg and H. G. Petersen, J. Chem. Phys. 91, 461 (1989)
    """

def get_mean(data):
    """Calculate mean
    """
    sum = 0
    for i in data:
        sum += i
    return sum*1.0/len(data)

def get_s2(data, ave):
    s2 = 0
    for i in data:
        s2 += (i-ave)*(i-ave)
    return s2*1.0

def split_data(data, seq):
    x = []
    if len(data)%seq == 0:
        res = len(data)/seq
    else:
        res = len(data)/seq + 1
    for i in range(res):
        if (i+1)*seq < len(data):
            x.append(data[i*seq: (i+1)*seq])
        else:
            x.append(data[i*seq:])
    return x

def block_average(fname, outfile):
    """\sqrt{\frac{c_{0}(n/m)}{n/m-1}}
    n : the size of sample
    m : the block size
    """

    x, y = parse_xvg(fname)
    ave = get_mean(y)
    print "Ave : ",ave,
    msd = math.sqrt(get_s2(y, ave)/(len(y)-1))
    print "Std Dev: ", msd
    o = open(outfile, 'w')
    for n in range(1, 20):
        ys = split_data(y, n)
        yn = []
        for i in ys:
            yn.append(get_mean(i))
        m = len(ys)
        if m > 1:
            c0 = math.sqrt(get_s2(yn, ave)/(m-1)/m)
        o.write("%15.4f%15.6f\n"%(n, c0))
    o.close()

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "ht:i:o:", \
                ["help", "type", "input", "output"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-t", "--type"):       
            t = arg
            if t == "xvg":
                func = parse_xvg
        elif opt in ("-o", "--output"):       
            outfile = arg
        elif opt in ("-i", "--input"):       
            infile = arg
    block_average(infile, outfile)
    
if __name__ == "__main__":
    main(sys.argv[1:]) 
