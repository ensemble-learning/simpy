"""
convert the simulation time to real time.
"""

import os
import argparse
import sys

def get_sim_time(args, n = 66):

    if args.mol:
        molname = args.mol[0]
    else:
        print "Exit! no molecular name given"
        sys.exit()

    if args.nmol:
        n = args.nmol[0]
    else:
        print "Exit! no molecular number given"
        sys.exit()

    flag = 0

    f = open("water.mol", "r")

    counter = 0
    for i in f :
        if "step" in i:
             step = int(i.strip().split()[0][4:])
             if counter == 0:
                 start = step
             counter += 1
        if molname in i:
            tokens = i.strip().split()
            nmol = int(tokens[0])
            if args.reactant:
                if nmol <= n:
                    flag = 1
                    break
            if args.product:
                if nmol >= n:
                    flag = 1
                    break

    f.close()

    if flag:
        return step, start
    else:
        print "Warning: not finish half reaction!"
        print "Only %d molecules have reacted"%nmol
        return 0, 0

def get_real_time(args, nstep, start):
    
    if args.fname:
        boostfile = args.fname[0] + ".bboost"
    
    if not os.path.exists(boostfile):
        realtime = nstep - start
    else:
        counter = 0
        realtime = 0.0
        f = open(boostfile, "r")
        for i in f:
            if counter > 0:
                tokens = i.strip().split()
                step = int(tokens[0])
                if step >= nstep:
                    break
                realtime += float(tokens[7])
            counter += 1
        f.close()
    return realtime

def get_half_time(args):
    nstep, start = get_sim_time(args)
    realtime = get_real_time(args,nstep, start)
    print nstep - start, realtime

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="water", nargs="?", help="simulation name")
    parser.add_argument("-reactant", action="store_true", help="input reactant")
    parser.add_argument("-product", action="store_true", help="input product")
    parser.add_argument("-mol", nargs=1, help="mol to analyze")
    parser.add_argument("-nmol", nargs=1, type=int,help="number of mols")
    #parser.add_argument("-vol", action="store_true", help="get the volume of the simulation box")
    args = parser.parse_args()
    get_half_time(args)

if __name__ == "__main__":
    main()
