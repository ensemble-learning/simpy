""" read the geo file and output to data (LAMMPS), geo and xyz file.
"""
from mytype import System, Molecule, Atom
from xyz import Xyz
from output_conf import toDump, toPdb
from block import xyzBlock

def movie():
    nframe = xyzBlock("movie.xyz", 242)

    print nframe

    if nframe > 1:
        for i in range(nframe+1):
            a = Xyz("output%05d.xyz"%i)
            b = a.parser()
            b.pbc = [14.547, 10.663, 14.384, 90.00, 91.06, 90.00] 
            b.step = i
            #toDump(b, "output%05d.dump"%i)
            toPdb(b, "output%05d.pdb"%i)
    """
    o = open("total.dump", "w")
    for i in range(nframe):
        f = open("output%05d.dump"%i, "r")
        for j in f:
            o.write(j)
        f.close()
    o.close()
    """
        

