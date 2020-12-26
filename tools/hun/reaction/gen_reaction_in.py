import os

# determine atom types
atomtypes = []
f = open('lammps.data', 'r')
for i in f:
    if i.strip().startswith('Masses'):
        break

for i in f:
    if i.strip().startswith('Atoms'):
        break
    else:
        tokens = i.strip().split()
        if len(tokens) > 0:
            atomtypes.append(tokens[-1])

o = open('lammps.in', 'w')
o.write('''units             real
atom_style        charge
boundary          p p p

read_data        lammps.data

reset_timestep   0

pair_style       reax/c control.reaxc safezone 1.8 mincap 180

#----Neighbor Section----#

neighbor        1.0 bin
neigh_modify    delay 0 every 10 check no

''')
o.write('pair_coeff      * * ffield')
for i in range(len(atomtypes)):
    o.write(' %s'%atomtypes[i])
o.write('\n')

o.write('''
fix             QEQ all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

compute         reax all pair reax/c

variable eb     equal c_reax[1]
variable ea     equal c_reax[2]
variable elp    equal c_reax[3]
variable emol   equal c_reax[4]
variable ev     equal c_reax[5]
variable epen   equal c_reax[6]
variable ecoa   equal c_reax[7]
variable ehb    equal c_reax[8]
variable et     equal c_reax[9]
variable eco    equal c_reax[10]
variable ew     equal c_reax[11]
variable ep     equal c_reax[12]
variable efi    equal c_reax[13]
variable eqeq   equal c_reax[14]

#--------Output info--------

thermo         1
thermo_style    custom step etotal ke pe temp press vol v_eb v_ea v_elp v_emol v_ev v_epen v_ecoa v_ehb v_et v_eco v_ew v_ep v_efi v_eqeq cella cellb cellc cellalpha cellbeta cellgamma pxx pyy pzz
thermo_modify   line multi

dump            100 all custom 1 dump.lammpstrj id type x y z vx vy vz
dump_modify     100 sort id

rerun           movie.xyz dump x y z box no format xyz

''')

o.close()

o = open('control.reaxc', 'w')
o.write('''simulation_name         lammps  ! output files will carry this name + their specific extension

tabulate_long_range     0 ! denotes the granularity of long range tabulation, 0 means no tabulation
energy_update_freq      0
remove_CoM_vel          500 ! remove the trans. & rot. vel around the CoM every 'this many' steps

nbrhood_cutoff          5.0  ! near neighbors cutoff for bond calculations in A
hbond_cutoff            7.5  ! cutoff distance for hydrogen bond interactions
bond_graph_cutoff       0.3  ! bond strength cutoff for bond graphs
thb_cutoff              0.001 ! cutoff value for three body interactions
q_err                   1e-6  ! average per atom error norm allowed in GMRES convergence

geo_format              0    ! 0: xyz, 1: pdb, 2: bgf
write_freq              1    ! write trajectory after so many steps
traj_compress           0    ! 0: no compression  1: uses zlib to compress trajectory output
traj_title              lammps  ! (no white spaces)
atom_info               1    ! 0: no atom info, 1: print basic atom info in the trajectory file
atom_forces             0    ! 0: basic atom format, 1: print force on each atom in the trajectory file
atom_velocities         0    ! 0: basic atom format, 1: print the velocity of each atom in the trajectory file
bond_info               1    ! 0: do not print bonds, 1: print bonds in the trajectory file
angle_info              0    ! 0: do not print angles, 1: print angles in the trajectory file 
''')
o.close()

o = open('molfrag.sh', 'w')
o.write('''
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
''')
o.close()

os.makedirs('data', exist_ok=True)

o = open('control.file', 'w')
o.write('''# Input files
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
''')
o.close()

o = open('Atom.dic', 'w')
o.write('''Li 6.941000
Na 22.990000
K 39.098000
H 1.007900
B 10.811000
C 12.011000
N 14.007000
O 15.999000
F 18.998000
P 30.974000
S 32.065000
Al 26.982000
Cl 35.453000
Ca 40.078000
Zn 65.390000
''')
o.close()

o = open('Cutoff.dic', 'w')
o.write('''Li    Li      1.4000
    Li    Na      1.4000
    Li     K      1.4000
    Li     H      1.4000
    Li     B      1.4000
    Li     C      1.4000
    Li     N      1.4000
    Li     O      1.4000
    Li     F      1.4000
    Li     P      1.4000
    Li     S      1.4000
    Li    Al      1.4000
    Li    Cl      1.4000
    Li    Ca      1.4000
    Li    Zn      1.4000
    Na    Na      1.4000
    Na     K      1.4000
    Na     H      1.4000
    Na     B      1.4000
    Na     C      1.4000
    Na     N      1.4000
    Na     O      1.4000
    Na     F      1.4000
    Na     P      1.4000
    Na     S      1.4000
    Na    Al      1.4000
    Na    Cl      1.4000
    Na    Ca      1.4000
    Na    Zn      1.4000
     K     K      1.4000
     K     H      1.4000
     K     B      1.4000
     K     C      1.4000
     K     N      1.4000
     K     O      1.4000
     K     F      1.4000
     K     P      1.4000
     K     S      1.4000
     K    Al      1.4000
     K    Cl      1.4000
     K    Ca      1.4000
     K    Zn      1.4000
     H     H      0.4000
     H     B      0.4000
     H     C      0.4000
     H     N      0.4000
     H     O      0.4000
     H     F      0.4000
     H     P      0.4000
     H     S      0.4000
     H    Al      0.4000
     H    Cl      0.4000
     H    Ca      0.4000
     H    Zn      0.4000
     B     B      0.4000
     B     C      0.4000
     B     N      0.4000
     B     O      0.4000
     B     F      0.4000
     B     P      0.4000
     B     S      0.4000
     B    Al      0.4000
     B    Cl      0.4000
     B    Ca      0.4000
     B    Zn      0.4000
     C     C      0.4000
     C     N      0.4000
     C     O      0.4000
     C     F      0.4000
     C     P      0.4000
     C     S      0.4000
     C    Al      0.4000
     C    Cl      0.4000
     C    Ca      0.4000
     C    Zn      0.4000
     N     N      0.4000
     N     O      0.4000
     N     F      0.4000
     N     P      0.4000
     N     S      0.4000
     N    Al      0.4000
     N    Cl      0.4000
     N    Ca      0.4000
     N    Zn      0.4000
     O     O      0.4000
     O     F      0.4000
     O     P      0.4000
     O     S      0.4000
     O    Al      0.4000
     O    Cl      0.4000
     O    Ca      0.4000
     O    Zn      0.4000
     F     F      0.4000
     F     P      0.4000
     F     S      0.4000
     F    Al      0.4000
     F    Cl      0.4000
     F    Ca      0.4000
     F    Zn      0.4000
     P     P      0.4000
     P     S      0.4000
     P    Al      0.4000
     P    Cl      0.4000
     P    Ca      0.4000
     P    Zn      0.4000
     S     S      0.4000
     S    Al      0.4000
     S    Cl      0.4000
     S    Ca      0.4000
     S    Zn      0.4000
    Al    Al      0.4000
    Al    Cl      0.4000
    Al    Ca      0.4000
    Al    Zn      0.4000
    Cl    Cl      0.4000
    Cl    Ca      0.4000
    Cl    Zn      0.4000
    Ca    Ca      0.4000
    Ca    Zn      0.4000
    Zn    Zn      0.4000
''')
o.close()
