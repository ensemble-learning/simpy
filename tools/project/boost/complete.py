"""
convert the simulation time to real time.
"""

import os
import argparse

def get_sim_time(n = 66):

    f = open("water.mol", "r")
    for i in f :
        if "step" in i:
             step = int(i.strip().split()[0][4:])
        if "H2O1" in i:
            tokens = i.strip().split()
            nwater = int(tokens[0])
            if nwater >= n:
                break
    f.close()
    
    return step

def get_real_time(nstep):
    
    if not os.path.exists("water.bboost"):
        realtime = nstep
    else:
        counter = 0
        realtime = 0.0
        f = open("water.bboost", "r")
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

def get_half_time():
    nstep = get_sim_time(33)
    realtime = get_real_time(nstep)
    print nstep, realtime

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="geo", nargs="?", help="geo file name")
    parser.add_argument("-c", action="store_true", help="convert the file to other formats (geo, xyz, gjf, lammps)")
    parser.add_argument("-pbc", action="store_true", help="using default pbc 5nm * 5nm * 5nm")
    parser.add_argument("-b", nargs=2, type=int, help="get the bond distance between a1, a2, a3")
    parser.add_argument("-a", nargs=3, type=int,help="get the angle of a1-a2-a3")
    parser.add_argument("-vol", action="store_true", help="get the volume of the simulation box")
    args = parser.parse_args()
    get_half_time()

if __name__ == "__main__":
    main()
