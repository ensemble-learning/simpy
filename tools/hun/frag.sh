# get xyz from XDATCAR
perl ~/soft/simpy/tools/vasp/xdat2xyz.pl
# get Lammps data from POSCAR
python ~/soft/simpy/lib/e_2_contcar.py
# copy data
cp lammps.data lammps
cp movie.xyz lammps
# cd folder
cd lammps
# run lammps
lammps -in lammps_input
# lmp to config.out
dump2config.py rerun.lmp ./anal/config.out
# lmp to reaxbonds.out
trj2reaxbonds.py lammps.trj ./anal/reaxbonds.out
# cd folder
cd anal
# if have H
python ~/soft/simpy/tools/hun/refine_bonds.py
# analysis the species
./molfrag.sh control.file


