import sys
import os
import socket
import argparse

LIB = ''

if socket.gethostname() == "cluster.hpc.org":
    LIB = "/home/chengtao/packages/simpy/simpy/lib"
elif socket.gethostname() == "tao-laptop":
    LIB = "/home/tao/Nutstore/code/simupy/lib"
elif socket.gethostname() == "atom.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "ion.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "fermion.wag.caltech.edu":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant12":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant1":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "giant3":
    LIB = "/net/hulk/home6/chengtao/soft/simpy/lib"
elif socket.gethostname() == "mu01":
    LIB = "/home/tcheng/soft/simpy/lib"

sys.path.insert(0 , LIB)

from cons import ELEMENT2MASS

def gen_atomdic(atoms):
    o = open("Atom.dic", "w")
    for i in atoms:
        o.write("%s %.6f\n"%(i, ELEMENT2MASS[i]))
    o.close()

def gen_cutoff(atoms):
    bo = 0.4
    o = open("Cutoff.dic", "w")
    for i in range(len(atoms)):
        for j in range(i, len(atoms)):
            o.write("%6s%6s%12.4f\n"%(atoms[i], atoms[j], bo))
    o.close()

def gen_molfrag():
    o = open("molfrag.sh", "w")
    o.write("""
# argument
if [ $# -ne 1 ]
then
  echo "Usage : molfrag.sh [control.file]"
  exit
fi

Logfile="log.molfrag"

echo "Molecular Fragment Analysis Toolkit"
echo "-----------------------------------"
echo ""
echo "Programmed by Hyungjun Kim (linus) at 2007"
echo ""

#echo "Molecular Fragment Analysis Toolkit" > $Logfile
#echo "-----------------------------------" >> $Logfile
#echo "" >> $Logfile
#echo "Programmed by Hyungjun Kim (linus) at 2007" >> $Logfile
#echo "" >> $Logfile

echo " *Creating Bond Change Logs."
#echo " *Creating Bond Change Logs." >> $Logfile
bondlog $1
#cp dataTATB/Bond.log dataTATB/Bond.log.bak
#cp dataTATB/Nbondchg.log dataTATB/Nbondchg.log.bak

echo ""
echo " *Reducing Bond Change Logs."
echo " (Deleting Redundant reactions)"

#echo "" >> $Logfile
#echo " *Reducing Bond Change Logs." >> $Logfile
#echo " (Deleting Redundant reactions)" >> $Logfile
netbondlog $1

echo ""
echo " *Molecular Fragment Analysis"

#echo "" >> $Logfile
#echo " *Molecular Fragment Analysis" >> $Logfile
molfrag $1
molstat $1

MolStatList=$(grep MolStatTotal $1)
MolStatList=${MolStatList#MolStatTotal}

fragtable $1 $MolStatList

echo ""
echo " *Rxn Analysis"

#echo "" >> $Logfile
#echo " *Molecular Fragment Analysis" >> $Logfile
rxnanal $1

echo "DONE"
#rm bondlog netbondlog molfrag rxnanal molstat fragtable
""")
    o.close()
    os.system("chmod +x molfrag.sh")

def gen_control():
    o = open("control.file", "w")
    o.write("""# Input files
# Type : Grasp, Reax, Reax2
Type            Grasp
Config          config.out
Bonds           reaxbonds.out

# Simulation Information (fs)
Timestep        0.5

# Dictionaries
Atoms           Atom.dic
Bond-CutOff     Cutoff.dic

# Analysis Parameters
# 0: Do not create MoleculeData file. 1: Create MoleculeData file.
MolDataFlag     1
# Time window for determining important bond breaking (ps)
TimeWindow      0.05

# Output files
MoleculeID              data/molid.out
MoleculeData    data/moldata.out

NofMolecules    data/nmol.out
MoleculeChanges data/molchg.out
MoleculeAtom    data/mol.atom
BondLog         data/Bond.log
NBondChange     data/Nbondchg.log
NetBonds        data/netBonds.out
RxnLog          data/rxn.log
TemperatureLog  data/temperature.log

# Output files for Molecular Statistics
MolStatPerStep  data/molstat.step
MolStatTotal    data/molstat.total
FragmentTable   data/fragtable

# Basic Settings for Memory Allocation (integer type)
MaxBondsPerAtom                         50
MaxAtomsPerMolecule                     2000
MaxMoleculeChangePerStep        100000
MaxMoleculeSpecies                      1000000
MaxMoleculeTrace                        10000

# Analyzing Frequency (integer type)
# Currently Frequency=1 case is only tested.
Frequency                       1
""")
    o.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-atoms", nargs="+", help="atom types")
    parser.add_argument("-d", nargs="+", help="dump file")
    parser.add_argument("-t", nargs="+", help="traj file")
    args = parser.parse_args()
    
    atoms = []
    if args.atoms:
        atoms = args.atoms
    else:
        exit("No atoms input!")
    
    if not os.path.exists("anal"):
        os.mkdir("anal")

    os.chdir("anal")
    if not os.path.exists("data"):
        os.mkdir("data")

    # generate the necessary input files
    gen_molfrag()
    gen_atomdic(atoms)
    gen_cutoff(atoms)
    gen_control()
    os.chdir("..")
    
    if 1:
        os.system("dump2config.py rerun.lmp ./anal/config.out")
        os.system("trj2reaxbonds.py lammps.trj ./anal/reaxbonds.out")
        os.chdir("anal")
        os.system("./molfrag.sh control.file")

if __name__ == "__main__":
    main()
