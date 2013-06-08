"""
lammps input templates
"""

MIN = """
#Lammps 2009 input file generate by DFF

units          real
atom_style     charge
boundary       p p p

read_data      lammps.data
#read_restart       add.rst

pair_style      reax/c control
#pair_style      reax/c NULL lgvdw yes

#----Neighbor Section----#
neighbor                1.0 bin
neigh_modify    delay 0 every 10 check no

#read_data      input.data
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

fix             QEQ all qeq/reax 1 0.0 10.0 1.0e-6 reax/c


thermo         1
thermo_style    custom step etotal ke pe temp press vol v_eb v_ea v_elp v_emol v_ev v_epen v_ecoa v_ehb v_et v_eco v_ew v_ep v_efi v_eqeq
thermo_modify   line multi
dump                    1 all custom 1 dump.lmp id type x y z vx vy vz
dump_modify             1 sort id

min_style       cg
minimize        0 1.0e-8 1000 1000
fix            9 all reax/c/bonds 1 1 1 bonds.reax
write_restart  min.rst
"""


NVT = """
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

fix             QEQ all qeq/reax 1 0.0 10.0 1.0e-6 reax/c

#--------Output info--------
thermo         400
thermo_style    custom step etotal ke pe temp press vol v_eb v_ea v_elp v_emol v_ev v_epen v_ecoa v_ehb v_et v_eco v_ew v_ep v_efi v_eqeq
thermo_modify   line multi
dump                    1 all custom 400 dump.lmp id type x y z vx vy vz
dump_modify             1 sort id

#--------Simulation--------
fix            101 all nvt temp 300 300 200       
fix            909 all reax/c/bonds 1 1 400 bonds.reax
timestep       0.25
run            80000
write_restart  nvt.rst

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


