""" read the geo file and output to data (LAMMPS), geo and xyz file.
"""
import sys, os
from mytype import System, Molecule, Atom
from g03 import G03Gjf
from output_conf import toXyz, toGeo, toGjf

usage = """gau2xyz g03logfiles"""
if len(sys.argv) < 2:
    print usage
else:
    for i in sys.argv[1:]:
        outfile = i.split(".")[0]
        if os.path.exists(i):
            a = G03Gjf(i)
            b = a.parser()
            b.geotag = "BIOGRF 200"
            toXyz(b, outfile+".xyz")
            toGeo(b, outfile+".geo")
            toGjf(b, outfile+".gjf.bak")
        else:
            print "missing file %s"%i
