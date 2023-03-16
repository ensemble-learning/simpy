#!/usr/bin/env python

import sys
import argparse
from ase import Atoms
from ase.units import Bohr, Hartree
from ase.io import read, write

SHE = 4.44
EV_2_HARTREE = 27.2114

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

def write_control(args):
    o = open('control', 'w')
    o.write(r"ion-species GBRV/$ID_pbe.uspp"+"\n")
    o.write(r"elec-ex-corr gga-PBE"+"\n")
    o.write(r"elec-cutoff 20"+"\n")
    #o.write(r"electronic-SCF energyDiffThreshold 3.6749e-7 nIterations 1000  # 1E-5 eV #default=1e-8 Ha"+"\n")
    o.write(r"elec-smearing Fermi 0.005"+"\n")
    o.write("\n")
    if args.spin:
        o.write("spintype z-spin" + "\n")
    o.write(r"symmetry-threshold 0.001"+"\n")
    o.write(r"coulomb-interaction Periodic"+"\n")
    o.write("\n")
    o.write(r"# input and output"+"\n")
    o.write(r"initial-state $VAR"+"\n")
    o.write(r"dump-name $VAR"+"\n")
    o.write(r"dump End State"+"\n")
    o.write(r"dump End Forces"+"\n")
    o.write(r"dump End Ecomponents"+"\n")
    o.write("\n")
    o.write(r"# solvation"+"\n")
    o.write(r"fluid LinearPCM"+"\n")
    o.write(r"pcm-variant CANDLE"+"\n")
    o.write(r"fluid-solvent H2O"+"\n")
    o.write(r"fluid-cation K+ 0.1"+"\n")
    o.write(r"fluid-anion F- 0.1"+"\n")
    o.write("\n")

    if args.cp:
        u = args.cp[0]
        target_u = -(SHE + u)/EV_2_HARTREE
        o.write(r"# constant potential" + "\n")
        o.write("target-mu %.10f"%target_u + "\n")

    if args.q:
        q = args.q[0]
        o.write("elec-initial-charge %.10f"%q+ "\n")

    o.close()

def write_input():
    o = open('jdft.in', 'w')
    o.write('include coords\n')
    o.write('include kpoints\n')
    o.write('include control\n')
    o.close()

def write_pbs():
    o = open('pbs', 'w')
    o.write("""#!/bin/bash
#PBS -N jdftx
#PBS -l nodes=1:ppn=32
#PBS -l walltime=90:00:00
#PBS -q batch
#PBS -S /bin/bash
#PBS -j oe

source /home/tcheng/intel/oneapi/setvars.sh intel64

nProcesses=2
JDFTx="/home/tcheng/soft/jdftx/jdftx-1.7.0/build/jdftx"
MPIRUN="mpirun --map-by node --bind-to none"
nAllocCores=$(wc -l < $PBS_NODEFILE)
nThreads=$(echo "$nAllocCores/$nProcesses" | bc)
hostname > $PBS_O_WORKDIR/$PBS_JOBNAME.machine
rm -rf /temp1/$USER/$PBS_JOBNAME
mkdir -p /temp1/$USER/$PBS_JOBNAME

cd $PBS_O_WORKDIR
$MPIRUN -n $nProcesses $JDFTx -c $nThreads -i jdft.in > $PBS_O_WORKDIR/jdft.out

""")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", default=["POSCAR"], nargs="*", help="inputfile")
    parser.add_argument("-cp", nargs=1, type=float,help="the voltage (in V)")
    parser.add_argument("-q", nargs=1, type=float,help="excess electrons compared to a neutral system")
    parser.add_argument("-spin", action='store_true', help="spin-polarized calculation")
    args = parser.parse_args()

    infile = args.infile[0]

    atoms = read(infile)
    write_coords(atoms)
    write_kpoints()
    write_control(args)
    write_input()
    write_pbs()

if __name__ == "__main__":
    main()
