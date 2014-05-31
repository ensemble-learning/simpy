"""
lammps input templates
"""

MIN = """
#Lammps 2009 input file generate by DFF

units          real
atom_style     charge
boundary       p p p

read_data      lammps.data

pair_style      %reax_potential%
pair_coeff      * * ffield %ffield_atoms%

#----Neighbor Section----#
neighbor                1.0 bin
neigh_modify    delay 0 every 10 check no


#----ReaxFF Energy Terms----#
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

fix             QEQ all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

thermo          1
thermo_style    custom step etotal ke pe temp press vol v_eb v_ea v_elp v_emol v_ev v_epen v_ecoa v_ehb v_et v_eco v_ew v_ep v_efi v_eqeq cella cellb cellc cellalpha cellbeta cellgamma
thermo_modify   line multi

#fix             201 all box/relax aniso 0.0 vmax 0.001

min_style       cg
minimize        0 1.0e-8 1000 1000
write_restart   min.rst

dump           100 all custom 1 dump.lammpstrj id type x y z vx vy vz
dump_modify    100 sort id
dump           101 all cfg 1 dump.*.cfg mass type xs ys zs vx vy vz fx fy fz
dump_modify    101 element %elements%

#fix            200 all reax/c/bonds 1 bonds.reax
#fix            201 all reax/c/species 1 1 1 species.out
run            1
"""

NVT = """
units             real
atom_style        charge
boundary          p p p

#read_data        lammps.data
read_restart      min.rst

reset_timestep   0

#pair_style      reax/c control
#pair_style      reax/c NULL lgvdw yes
pair_style       %reax_potential%

#----Neighbor Section----#

neighbor        1.0 bin
neigh_modify    delay 0 every 10 check no

pair_coeff      * * ffield %ffield_atoms%
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

thermo         400
thermo_style    custom step etotal ke pe temp press vol v_eb v_ea v_elp v_emol v_ev v_epen v_ecoa v_ehb v_et v_eco v_ew v_ep v_efi v_eqeq cella cellb cellc cellalpha cellbeta cellgamma
thermo_modify   line multi

dump            100 all custom 400 dump.lammpstrj id type x y z vx vy vz
dump_modify     100 sort id

#--------Analysis-----------

fix             200 all reax/c/bonds 400 bonds.reax
#fix             201 all reax/c/species 1 1 400 species.out

#--------Simulation--------

velocity        all create 300.0 4928459 rot yes dist gaussian
fix             401 all nvt temp 300 300 50       
timestep        0.25
run             80000
write_restart   nvt.rst

"""


NPT = """
units          real
atom_style     charge
boundary       p p p

read_data      lammps.data
#read_restart       add.rst

#pair_style      reax/c control
pair_style      reax/c NULL lgvdw yes

#----Neighbor Section----#
neighbor                1.0 bin
neigh_modify    delay 0 every 10 check no

pair_coeff      * * ffield %ffield_atoms%

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

variable       density equal mass(all)*10/(6.0221417930*vol)

fix             QEQ all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

#--------Output info--------
thermo         400
thermo_style    custom step etotal ke pe temp press vol v_density v_eb v_ea v_elp v_emol v_ev v_epen v_ecoa v_ehb v_et v_eco v_ew v_ep v_efi v_eqeq
thermo_modify   line multi
dump                    1 all custom 400 dump.lmp id type x y z vx vy vz
dump_modify             1 sort id

#--------Simulation--------
fix            101 all npt temp 300 300 100 iso 1.0 1.0 25       
fix            909 all reax/c/bonds 1 1 400 bonds.reax
timestep       0.25
run            80000
write_restart  npt.rst

"""

CONTROL = """simulation_name         lammps  ! output files will carry this name + their specific extension

tabulate_long_range     0 ! denotes the granularity of long range tabulation, 0 means no tabulation
energy_update_freq      0
remove_CoM_vel          500 ! remove the trans. & rot. vel around the CoM every 'this many' steps

nbrhood_cutoff          5.0  ! near neighbors cutoff for bond calculations in A
hbond_cutoff            3.5  ! cutoff distance for hydrogen bond interactions
bond_graph_cutoff       0.3  ! bond strength cutoff for bond graphs
thb_cutoff              0.001 ! cutoff value for three body interactions
q_err                   1e-6  ! average per atom error norm allowed in GMRES convergence

geo_format              0    ! 0: xyz, 1: pdb, 2: bgf
write_freq              800    ! write trajectory after so many steps
traj_compress           0    ! 0: no compression  1: uses zlib to compress trajectory output
traj_title              lammps  ! (no white spaces)
atom_info               1    ! 0: no atom info, 1: print basic atom info in the trajectory file
atom_forces             0    ! 0: basic atom format, 1: print force on each atom in the trajectory file
atom_velocities         0    ! 0: basic atom format, 1: print the velocity of each atom in the trajectory file
bond_info               1    ! 0: do not print bonds, 1: print bonds in the trajectory file
angle_info              0    ! 0: do not print angles, 1: print angles in the trajectory file 
"""
