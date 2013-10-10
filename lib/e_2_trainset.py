
import os
import argparse
#from utilities import parseBlock
from mytype import System, Molecule, Atom
from trainset import Trainset
#from output_conf import toGeo, toXyz, toReaxLammps, toPdb, toGjf

def large_than_10():
    a = Trainset()
    for i in range(len(a.energy_lines)):
        md = a.energy[0][i]
        qm = a.energy[1][i]
        if abs(md - qm) > 10.0:
            print a.energy_lines[i].strip()
        else:
            print abs(md - qm)

if __name__ == "__main__":
    large_than_10()



