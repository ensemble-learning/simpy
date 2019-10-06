#!/usr/bin/env python

import sys
import argparse
from ase import Atoms
from ase.units import Bohr, Hartree
from ase.io import read, write

SHE = 4.44
EV_2_HARTREE = 27.2114

def read_atoms():
    if len(sys.argv) > 1: 
        atoms = read(sys.argv[1])
    else:
        atoms = read("POSCAR")
    return atoms

def write_coords(atoms):
    # Add lattice info
    R = atoms.get_cell()/Bohr
    o = open('coords', 'w')
    o.write('lattice \\\n')
    for i in range(3):
        for j in range(3):
            o.write('%f '%(R[j,i]))
        if(i!=2):
            o.write('\\')
        o.write('\n')
    # Add ion info
    atomPos = [x/Bohr for x in list(atoms.get_positions())]
    atomNames = atoms.get_chemical_symbols()
    o.write('\ncoords-type cartesian\n')
    fix_flag = 1 # 1: movable; 0: fixed
    for i in range(len(atomPos)):
        o.write('ion %s %f %f %f \t %d\n'%(atomNames[i], atomPos[i][0], atomPos[i][1], atomPos[i][2], fix_flag))
    del i
    # Construct most of the input file
    o.write('\n')
    o.close()

def write_kpoints():
    o = open('kpoints', 'w')
    o.write('kpoint-folding 1 1 1\n')
    o.close()

def write_control():
    o = open('control', 'w')
    o.write(r"""ion-species GBRV/$ID_pbe.uspp
elec-ex-corr gga-PBE
elec-cutoff 20 
electronic-SCF energyDiffThreshold 3.6749e-7   # 1E-5 eV #default=1e-8 Ha
elec-smearing Fermi 0.005

symmetry-threshold 0.001
coulomb-interaction Periodic

# input and output
initial-state $VAR
dump-name $VAR
dump End State
dump End Forces
dump End Ecomponents

# solvation
fluid LinearPCM
pcm-variant CANDLE
fluid-solvent H2O
fluid-cation K+ 0.1
fluid-anion F- 0.1

""")
    o.close()

def write_control_const_u(args):
    u = args.cp[0]
    target_u = -(SHE + u)/EV_2_HARTREE
    
    o = open('control', 'w')
    o.write(r"""ion-species GBRV/$ID_pbe.uspp
elec-ex-corr gga-PBE
elec-cutoff 20 
electronic-SCF energyDiffThreshold 3.6749e-7   # 1E-5 eV #default=1e-8 Ha
elec-smearing Fermi 0.005

symmetry-threshold 0.001
coulomb-interaction Periodic

# input and output
initial-state $VAR
dump-name $VAR
dump End State
dump End Forces
dump End Ecomponents

# solvation
fluid LinearPCM
pcm-variant CANDLE
fluid-solvent H2O
fluid-cation K+ 0.1
fluid-anion F- 0.1

# constant potential
target-mu %.10f 

"""%target_u)
    o.close()

def write_input():
    o = open('jdft.in', 'w')
    o.write('include coords\n')
    o.write('include kpoints\n')
    o.write('include control\n')
    o.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", default="POSCAR", nargs="*", help="inputfile")
    parser.add_argument("-cp", nargs=1, type=float,help="the voltage (in V)")
    args = parser.parse_args()

    infile = args.infile
    atoms = read(infile)
    write_coords(atoms)
    write_kpoints()
    if args.cp:
        write_control_const_u(args)
    else:
        write_control()
    write_input()

if __name__ == "__main__":
    main()
