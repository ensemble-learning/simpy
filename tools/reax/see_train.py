"""
Format and compare the fit error
"""

import os
import argparse

def parse_train(fname):
    data = []
    f = open(fname, "r")

    flag = 1
    while (flag):
        flag = 0
        for i in f:
            flag = 1
            if "GEOMETRY" in i:
                sec = "GEOMETRY"
                break
            elif "ENERGY" in i:
                sec = "ENERGY"
                break
            elif "CHARGES" in i:
                sec = "CHARGES"
                break
            else:
                pass

        for i in f:
            flag = 1
            if "Weight" in i:
                pass
            elif "Total Mean" in i:
                break
            else:
                if sec == "ENERGY":
                    tokens = i.strip().split()
                    qm = float(tokens[-5])
                    ff = float(tokens[-4])
                    err = float(tokens[-3])
                    wt = float(tokens[0])
                    reac = " ".join(tokens[2:-5])
                    data.append([wt, qm, ff, err, reac])
                elif sec == "CHARGES":
                    tokens = i.strip().split()
                    qm = float(tokens[2])
                    ff = float(tokens[3])
                    err = float(tokens[4])
                    #wt = float(tokens[0])
                    wt = 1.0
                    reac = tokens[1]
                    data.append([wt, qm, ff, err, reac])
    return data
        
def output(data):
    for i in data:
        print "%5.2f"%i[0],
        print "%40s"%i[4],
        print "%8.2f"%i[1],
        print "%8.2f"%i[2],
        print "%8.2f"%i[3]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="trainset.err.best", nargs="?", help="trainset.err file")
    parser.add_argument("-c", action="store_true", help="convert the file to other formats (geo, xyz, gjf, lammps)")
    #parser.add_argument("-pbc", action="store_true", help="using default pbc 5nm * 5nm * 5nm")
    #parser.add_argument("-b", nargs=2, type=int, help="get the bond distance between a1, a2, a3")
    #parser.add_argument("-a", nargs=3, type=int,help="get the angle of a1-a2-a3")
    #parser.add_argument("-vol", action="store_true", help="get the volume of the simulation box")

    args = parser.parse_args()

    if args.c:
        fname = "trainset.err.initial"
        data0 = parse_train(fname)
        fname = "trainset.err.best"
        data1 = parse_train(fname)
        assert len(data0) == len(data1)
        for i in range(len(data0)):
            data1[i][3] = data0[i][2]
        output(data1)
    else:
        fname = args.fname
        assert os.path.exists(fname)
        data = parse_train(fname)
        output(data)
    
if __name__ == "__main__":
    main()

