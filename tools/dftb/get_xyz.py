"""get x, y and z coordinations from geo_end.xyz
"""
import os
import sys

def get_xyz():
    if not os.path.exists("geo_end.xyz"):
        print "Error: No geo_end.xyz file in current path!"
        sys.exit()
    if os.path.exists("lammps.xyz"):
        print "Warning: Overwrite lammps.xyz !"

    f = open("geo_end.xyz", "r")
    o = open("lammps.xyz", "w")

    for i in f:
        tokens = i.strip().split()
        #@note: not a strick way to get the coords
        if len(tokens) == 8:
            line = " ".join(tokens[:4]) + "\n"
        else:
            line = i
        o.write(line)
    o.close()
    f.close()

if __name__ == "__main__":
    get_xyz()
