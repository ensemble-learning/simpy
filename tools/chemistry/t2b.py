#!/usr/bin/env python

"""
usage: t2b.py tafel_slope
convert the tafel slope to apparent transfer coefficient
ref: acscatal.7b03142
"""

import sys

# a = kbTln(10)/e mV/dec
a = 59

if len(sys.argv) > 1:
    b = float(sys.argv[1])
    atc = a/b # apparent transfer coefficient
    print("apparent transfer coefficient (atc) = %.3f"%atc)
    if atc < 1:
        print("atc at TS = %.3f after 0 electron transfer"%atc)
    if atc >= 1 and atc < 2:
        print("atc at TS = %.3f after 1 electron transfer"%(atc-1))
    if atc >= 2 and atc < 3:
        print("atc at TS = %.3f after 2 electron transfer"%(atc-2))
    if atc >= 3 and atc < 4:
        print("atc at TS = %.3f after 2 electron transfer"%(atc-3))
    if atc >= 4 and atc < 5:
        print("atc at TS = %.3f after 2 electron transfer"%(atc-4))
    if atc >= 5 and atc < 6:
        print("atc at TS = %.3f after 2 electron transfer"%(atc-5))
    if atc >= 6 and atc < 7:
        print("atc at TS = %.3f after 2 electron transfer"%(atc-6))
    if atc >= 7 and atc < 8:
        print("atc at TS = %.3f after 2 electron transfer"%(atc-7))
else:
    print(__doc__)

