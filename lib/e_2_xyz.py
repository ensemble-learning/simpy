""" read the geo file and output to data (LAMMPS), geo and xyz file.
"""
from mytype import System, Molecule, Atom
from g03 import G03LogConf
from output_conf import toXyz

for i in range(21):
    testfile = "scan%02d.log"%i
    a = G03LogConf(testfile)
    b = a.parser()
    toXyz(b, "scan%02d.xyz"%i)
