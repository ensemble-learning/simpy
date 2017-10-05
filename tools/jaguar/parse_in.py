"""
Parse the Jaguar in file
"""
import sys, os
import shutil
import re

class System:
    def __init__(self,):
        self.orbitals = []
    def get_homo(self,):
        for i in self.orbitals:
            if i.occupy > 0.0:
                id = i.orbital_id
            else:
                break
        for i in range(id-9, id+10):
            if i < id:
                print "Homo-%d"%(id-i),
            elif i == id:
                print "Homo  ",
            elif i == id + 1:
                print "Lumo  ",
            else:
                print "Lumo-%d"%(i-id-1),
            print i, 
            print self.orbitals[i-1].occupy, 
            print "%12.6f"%self.orbitals[i-1].energy,
            print "%12.6f"%((self.orbitals[i-1].energy - self.orbitals[id-1].energy)*627.509)

class OrbitalEnergy:
    def __init__(self,):
        self.orbital_id = 0
        self.occupy = 0.0
        self.energy = 0.0

def parse_in(fname, system):
    f = open(fname, "r")
    for i in f:
        if i.strip().startswith("&guess"):
            break
    for i in f:
        if i.strip().startswith("&"):
            break
        else:
            if "Orbital Energy" in i:
                oe = OrbitalEnergy()
                tokens = i.strip().split()
                oe.orbital_id = int(tokens[0])
                oe.energy = float(tokens[3])
                oe.occupy = float(tokens[5])
                system.orbitals.append(oe)
    
def usage():
    sys.stdout.write("python parse_in.py infile\n")

def main():
    if len(sys.argv) < 2:
        usage()
    else:
        system = System()
        fname = sys.argv[1]
        parse_in(fname, system)
        system.get_homo()

if __name__ == "__main__":
    main()
