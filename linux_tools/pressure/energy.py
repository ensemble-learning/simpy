import os
for j in range(1,21):
    os.system("mkdir c%d"%j)
    os.system("cp c%d.gro c%d"%(j,j))
    shfile = open("c%d/pressure.sh"%j ,'w')
    shfile.write("""#!/bin/bash
#PBS -l nodes=1:ppn=8
cd $PBS_O_WORKDIR
""")
    for i in range(31):
        cutoff = 10 + i
        o = open("c%d/a%d.mdp"%(j,cutoff), 'w')
        o.write("""title           = viscosity for spc 1728


cpp             = /usr/bin/cpp
integrator      = md
dt              = 0.0010 ; ps !
nsteps          = 0 ;
nstcomm         = 1
nstxout         = 0 ; collect data every 1.0 ps
nstvout         = 0
nstfout         = 0
nstlist         = 5
ns_type         = grid
rlist           = %4.1f
coulombtype     = cut-off 
rcoulomb        = %4.1f
rvdw            = %4.1f
DispCorr        = EnerPres
fourierspacing  = 0.12
fourier_nx      = 0
fourier_ny      = 0
fourier_nz      = 0
pme_order       = 4
ewald_rtol      = 1e-5
optimize_fft    = yes
; Berendsen temperature coupling is on in two groups
Tcoupl          = berendsen
tc_grps         = System
tau_t           = 0.1
ref_t           = 273.0
; Pressure coupling is on
Pcoupl          = no ;berendsen ;parrinello-rahman  
pcoupltype      = isotropic
tau_p           = 0.5
compressibility = 4.5e-5
ref_p           = 1.000
; Generate velocites is on at 300 K.
gen_vel         = no ;yes 
gen_temp        = 300.0
gen_seed        = 94823
; Constrains of Bonds
constraints      = all-bonds
constraint_algorithm = shake
shake_tol        = 0.0001
; Viscosity calculations
;cos_acceleration = 0.06 
"""%(float(cutoff) / 10, float(cutoff) / 10, float(cutoff) / 10))
        o.close()
        shfile.write("grompp_mpi -np 8 -f a%d.mdp -c c%d.gro -p *.top -o equil\n"%(cutoff, j))
        shfile.write("/opt/openmpi/bin/mpirun -np 8 mdrun_mpi -s equil -g c%da%d.log\n"%(j,cutoff))
    shfile.close()
