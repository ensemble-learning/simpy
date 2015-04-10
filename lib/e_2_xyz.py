"""Read the geo file and output to data (LAMMPS), geo and xyz file.
@todo: 
2015-04-06: Print doc string if no arguments supplied.
"""
import sys, os
import stat
import time
import argparse
import subprocess
import ConfigParser
from mytype import System, Molecule, Atom
from xyz import Xyz
from output_conf import toDump, toPdb
from block import xyzBlock

class XyzInfo():
    def __init__(self,):
        self.natoms = 0
        self.pbc = [0.0, 0.0, 0.0, 90.0, 90.0, 90.0]
        self.restart = 0

def ConfigSectionMap(config, section):
    dict1 = {}
    options = config.options(section)
    for option in options:
        try:
            dict1[option] = config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            sys.stdout.write("exception on %s!\n" % option)
            dict1[option] = None
    return dict1

def replicate(xyz_info, args):
    fname = args.f[0]
    nlines = xyz_info.natoms + 2
    nframe = xyzBlock(fname, nlines)

    na, nb, nc = args.nbox
    if nframe > 1:
        for i in range(nframe+1):
            if not xyz_info.restart:
                a = Xyz("output%05d.xyz"%i)
                b = a.parser()
                b.pbc = xyz_info.pbc
                b.step = i
                toPdb(b, "output%05d.pdb"%i)
            # replicate the box using external code (genconf in Gromacs)
            subprocess.Popen(["genconf", "-f", "output%05d.pdb"%i, 
                              "-o", "s_%05d.pdb"%i, "-nbox",
                                "%d"%na, "%d"%nb, "%d"%nc],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    subprocess.Popen(["cat s_*.pdb > movie_%d_%d_%d.pdb"%(na, nb, nc)],
                     shell=True)
    # clear up
    o = open("clear.sh", "w")
    o.write("rm output*.xyz\n")
    o.write("rm output*.pdb\n")
    o.write("rm s_*.pdb\n")
    o.close()
    st = os.stat('clear.sh')
    os.chmod('clear.sh', st.st_mode | stat.S_IEXEC)

def read_info(xyz_info, config):
    tokens = ConfigSectionMap(config, "XYZ")['pbc']
    tmp = [float(i) for i in tokens.strip().split()]
    for i in range(len(tmp)):
        xyz_info.pbc[i] = tmp[i]
    xyz_info.natoms = int(ConfigSectionMap(config, "XYZ")['natoms'])
    flag = ConfigSectionMap(config, "XYZ")['restart']
    flag = flag.capitalize()
    if flag == "YES":
        xyz_info.restart = 1
    else:
        xyz_info.restart = 0

def write_config():
    """ Write an example ini file.
    """
    o = open("test.ini", "w")
    o.write("""[XYZ]
pbc = 11.155 11.155 50.0 90.0 90.0 120.0
natoms = 194
restart = yes
""")
    o.close()

def main(args):
    xyz_info = XyzInfo()

    if args.f:
        fname = args.f[0]
    else:
        fname = "movie.xyz"
        sys.stderr.write("Using default xyz file: movie.xyz\n")

    if args.s:
        config = ConfigParser.ConfigParser() 
        config.read(args.s[0])
        read_info(xyz_info, config)

    if args.nbox:
        if not args.s:
            sys.stderr.write("No configure file found!\n")
            sys.stderr.write("Prepare a sample configure file (test.ini)\n")
            write_config()
            sys.exit()
        else:
            replicate(xyz_info, args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", nargs=1, help="xyz file name")
    parser.add_argument("-s", nargs=1, help="configure file")
    parser.add_argument("-nbox", nargs=3, type=int, help="replicate the box in a, b and c")
    args = parser.parse_args()

    main(args)
    
