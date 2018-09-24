#!/usr/bin/env python

"""
@note: require MCCCS Towhee 7-2-0
@usage: pdb2towheecoords pdbfile nmol
"""

import sys, os, socket
LIB = ''

print socket.gethostname()

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simpy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simupy/lib"
elif socket.gethostname() == "atom.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "ion.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant12":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif "comet" in socket.gethostname():
    LIB = "/home/tcheng/soft/simpy/lib"
elif socket.gethostname() == "tao-Precision-Tower-3420":
    LIB = "/home/tao/Soft/simpy/lib"
elif socket.gethostname() == "zwicky":
    LIB = "/home/tcheng/Soft/simpy/lib"
elif "onyx" in socket.gethostname():
    LIB = "/p/home/taocheng/src/simpy/lib"
elif socket.gethostname() == "tao-ThinkCentre-M79":
    LIB = "/home/tao/Soft/simpy/lib"
elif socket.gethostname() == "tao-Precision-Tower-3420-ubuntu":
    LIB = "/home/tao/soft/simpy/lib"
elif "stampede2" in socket.gethostname():
    LIB = "/home1/04076/tg833760/soft/simpy/lib"

sys.path.insert(0 , LIB)

from mytype import System, Molecule, Atom
from pdb import Pdb
from output_conf import toTowheecoords
from utilities import lattice2v
from cons import ELEMENT2MASS, ELEMENT2ATN, A2Bohr

UFF_path = "/home/tao/data/soft/towhee/towhee-7.2.0/ForceFields/towhee_ff_UFF"

UFF ={"Pt":"Pt4+2", "Cl":"Cl", "O":"O_3", 
      "H":"H_", "C":"C_3", "N":"N_3", "Fe":"Fe3+2",
      "Cu":"Cu3+1", "Au":"Au4+3"}

water = """#water
input_style
'basic connectivity map'
nunit
3
nmaxcbmc
3
lpdbnames
F
forcefield
'UFF'
charge_assignment
'manual'
unit ntype
1    'H_'  0.4d0
vibration
1
2
improper
0
unit ntype
2    'O_3' -0.8d0
vibration
2
1 3
improper
0
unit ntype
3    'H_'  0.4d0
vibration
1
2
improper
0
"""

if len(sys.argv) > 2:
    nmol = int(sys.argv[2])
    a = Pdb(sys.argv[1])
    b = a.parser()
    b.assignAtomTypes()
    b.assignEleTypes()
    toTowheecoords(b)

    xx, xy, xz, yy, yz, zz = lattice2v(b.pbc)
    la = [xx, 0.0, 0.0]
    lb = [xy, yy, 0.0]
    lc = [xz, yz, zz]
    print la
    print lb
    print lc
    
    o = open("towhee_input", "w")
    o.write("inputformat\n")
    o.write("'Towhee'\n")
    o.write("ensemble\n")
    o.write("'uvt'\n")
    o.write("temperature\n")
    o.write("256.0d0\n")
    o.write("nmolty\n")
    o.write("2\n")
    o.write("nmolectyp\n")
    o.write("1 %d\n"%nmol)
    o.write("chempot\n")
    o.write("#10422.0d0\n")
    o.write("-14359.426d0\n")
    o.write("143509.426d0\n")
    o.write("numboxes\n")
    o.write("1\n")
    o.write("stepstyle\n")
    o.write("'moves'\n")
    o.write("nstep\n")
    o.write("%d\n"%(nmol*20))
    o.write("printfreq\n")
    o.write("10   \n")
    o.write("blocksize\n")
    o.write("20\n")
    o.write("moviefreq\n")
    o.write("100000\n")
    o.write("backupfreq\n")
    o.write("1000\n")
    o.write("runoutput\n")
    o.write("'full'\n")
    o.write("pdb_output_freq\n")
    o.write("100000\n")
    o.write("louthist\n")
    o.write("T\n")
    o.write("hist_label\n")
    o.write(" 1\n")
    o.write("hist_suffix\n")
    o.write("'d'\n")
    o.write("hist_nequil\n")
    o.write(" 10\n")
    o.write("histcalcfreq\n")
    o.write("5\n")
    o.write("histdumpfreq\n")
    o.write("10\n")
    o.write("pressure_virial_freq\n")
    o.write("20\n")
    o.write("trmaxdispfreq\n")
    o.write("1000\n")
    o.write("volmaxdispfreq\n")
    o.write("1000\n")
    o.write("potentialstyle\n")
    o.write("'internal'\n")
    o.write("ffnumber\n")
    o.write("1\n")
    o.write("ff_filename\n")
    o.write("%s\n"%UFF_path)
    o.write("classical_potential\n")
    o.write("'UFF 12-6'       \n")
    o.write("classical_mixrule\n")
    o.write("'Geometric'\n")
    o.write("lshift\n")
    o.write(".false.\n")
    o.write("ltailc\n")
    o.write(".true.\n")
    o.write("rmin  \n")
    o.write("0.6d0 \n")
    o.write("rcut  \n")
    o.write("4.0d0\n")
    o.write("rcutin \n")
    o.write("4.0d0 \n")
    o.write("electrostatic_form\n")
    o.write("'coulomb'\n")
    o.write("coulombstyle\n")
    o.write("'ewald_fixed_kmax'\n")
    o.write("kalp \n")
    o.write("5.6d0\n")
    o.write("kmax\n")
    o.write("5\n")
    o.write("dielect\n")
    o.write("1.0d0\n")
    o.write("linit   \n")
    o.write("T\n")
    o.write("initboxtype\n")
    o.write("'dimensions'\n")
    o.write("initstyle\n")
    o.write("'coords'\n")
    o.write("'full cbmc'\n")
    o.write("initlattice\n")
    o.write("'none' \n")
    o.write("'none' \n")
    o.write("initmol\n")
    o.write("1\n")
    o.write("0\n")
    o.write("inix iniy iniz\n")
    o.write("3    3    9   \n")
    o.write("hmatrix\n")
    for ii in la:
        o.write("%.6fd0 "%ii)
    o.write("\n")
    for ii in lb:
        o.write("%.6fd0 "%ii)
    o.write("\n")
    for ii in lc:
        o.write("%.6fd0 "%ii)
    o.write("\n")
    o.write("pmuvtcbswap\n")
    o.write("0.25d0\n")
    o.write("          pmuvtcbmt\n")
    o.write("          0.0d0 1.0d0\n")
    o.write("pm1boxcbswap\n")
    o.write("0.0d0\n")
    o.write("          pm1cbswmt\n")
    o.write("          0.0d0 1.0d0\n")
    o.write("pmavb1\n")
    o.write("0.0d0\n")
    o.write("          pmavb1in\n")
    o.write("          0.5d0\n")
    o.write("          pmavb1mt\n")
    o.write("          1.0d0 1.0d0\n")
    o.write("          pmavb1ct\n")
    o.write("          1.0d0 1.0d0\n")
    o.write("          1.0d0 1.0d0\n")
    o.write("          avb1rad\n")
    o.write("          4.0d0\n")
    o.write("pmavb2\n")
    o.write("0.0d0\n")
    o.write("          pmavb2in\n")
    o.write("          0.5\n")
    o.write("          pmavb2mt\n")
    o.write("          1.0d0 1.0d0\n")
    o.write("          pmavb2ct\n")
    o.write("          1.0d0 1.0d0\n")
    o.write("          1.0d0 1.0d0\n")
    o.write("          avb2rad\n")
    o.write("          4.0d0\n")
    o.write("pmavb3\n")
    o.write("0.0d0\n")
    o.write("          pmavb3mt\n")
    o.write("          1.0d0 1.0d0\n")
    o.write("          pmavb3ct\n")
    o.write("          1.0d0 1.0d0\n")
    o.write("          1.0d0 1.0d0\n")
    o.write("          avb3rad\n")
    o.write("          4.0d0\n")
    o.write("pmcb\n")
    o.write("0.50d0\n")
    o.write("          pmcbmt\n")
    o.write("          0.0d0 1.0d0\n")
    o.write("          pmall\n")
    o.write("          0.0d0 0.0d0\n")
    o.write("pmback\n")
    o.write("0.0d0\n")
    o.write("          pmbkmt\n")
    o.write("          0.0d0 1.0d0\n")
    o.write("pmpivot\n")
    o.write("0.0d0\n")
    o.write("          pmpivmt\n")
    o.write("          0.0d0 1.0d0\n")
    o.write("pmconrot\n")
    o.write("0.0d0\n")
    o.write("          pmcrmt\n")
    o.write("          0.0d0 1.0d0\n")
    o.write("pmcrback\n")
    o.write("0.0d0\n")
    o.write("          pmcrbmt\n")
    o.write("          0.0d0 1.0d0\n")
    o.write("pmplane\n")
    o.write("0.0d0\n")
    o.write("          pmplanebox\n")
    o.write("          1.0d0\n")
    o.write("          planewidth\n")
    o.write("          3.0d0\n")
    o.write("pmrow\n")
    o.write("0.0d0\n")
    o.write("          pmrowbox\n")
    o.write("          1.0d0\n")
    o.write("          rowwidth\n")
    o.write("          3.0d0\n")
    o.write("pmtraat\n")
    o.write("0.0d0\n")
    o.write("          pmtamt\n")
    o.write("          1.0d0 0.0d0\n")
    o.write("          rmtraa\n")
    o.write("          0.5d0\n")
    o.write("          tatraa\n")
    o.write("          0.5d0\n")
    o.write("pmtracm\n")
    o.write("0.75d0\n")
    o.write("          pmtcmt\n")
    o.write("          0.0d0 1.0d0\n")
    o.write("          rmtrac\n")
    o.write("          0.5d0\n")
    o.write("          tatrac\n")
    o.write("          0.5d0\n")
    o.write("pmrotate\n")
    o.write("1.0d0\n")
    o.write("          pmromt\n")
    o.write("          0.0d0 1.0d0\n")
    o.write("          rmrot\n")
    o.write("          0.05d0\n")
    o.write("          tarot\n")
    o.write("          0.5d0\n")
    o.write("cbmc_formulation\n")
    o.write("'Martin and Siepmann 1999 + Martin and Thompson 2004'\n")
    o.write("cbmc_setting_style\n")
    o.write("'default ideal'\n")
    o.write("#substract\n")
    o.write("input_style\n")
    o.write("'basic connectivity map'\n")
    o.write("nunit\n")
    n = 0
    eles = []
    f = open("towhee_coords", "r")
    for i in f:
        tokens = i.strip().split()
        if len(tokens) > 3:
            ele = tokens[-1]
            eles.append(ele)
            n += 1
    f.close()
    o.write("%d\n"%n)
    o.write("nmaxcbmc\n")
    o.write("%d\n"%n)
    o.write("lpdbnames\n")
    o.write("F\n")
    o.write("forcefield\n")
    o.write("'UFF'\n")
    o.write("charge_assignment\n")
    o.write("'manual'\n")
    for i in range(len(eles)):
        o.write("unit ntype qqatom\n")
        o.write("%d    '%s'   0.0d0\n"%(i+1, UFF[eles[i]]))
        o.write("vibration\n")
        o.write("0\n")
        o.write("improper\n")
        o.write("0\n")
    o.write(water)
    #os.system("towhee")

else:
    pass
