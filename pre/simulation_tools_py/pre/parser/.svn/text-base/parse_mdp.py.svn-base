import os

MDP = r"""cpp             = /usr/bin/cpp
integrator      = md
dt              = 0.0020 ; ps !
nsteps          = %STEPS% ;
PBC             = no
nstcomm         = 1
nstxout         = 1000 ; collect data every 1.0 ps
nstvout         = 0
nstfout         = 0
nstlist         = 5
ns_type         = simple
comm_mode       = ANGULAR
rlist           = 3
coulombtype     = cut-off
rcoulomb        = 3
rvdw            = 3
DispCorr        = no
fourierspacing  = 0.12
fourier_nx      = 0
fourier_ny      = 0
fourier_nz      = 0
pme_order       = 4
ewald_rtol      = 1e-5
optimize_fft    = yes
; Berendsen temperature coupling is on in two groups
Tcoupl          = nose-hoover
tc_grps         = System
tau_t           = 0.1
ref_t           = %TEMP%
; Pressure coupling is on
Pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 0.5
compressibility = 4.5e-5
ref_p           = 1.000
; Generate velocites is on at 300 K.
gen_vel         = no
gen_temp        = 298.0
gen_seed        = 94823
; Constrains of Bonds
constraints      = hbonds
constraint_algorithm = lincs
shake_tol        = 0.0001
"""

TEMP = {
"01_2-propanol"  : 298.0,
"02_benzene" :  298.0,
"03_butane" :  273.0,
"0x_cyclohexane" :  298.0,
"04_ethane"  : 185.0,
"05_ethanol"  :  298.0,
"06_isobutane" :  298.0,
"07_methanol" :  298.0,
"08_phenol" :  298.0,
"09_propane"  :  230.0,
"10_propanol" :  298.0,
"11_propene" :   225.0,
"12_t-BuOH"  : 298.0,
"13_trans-2-butene" :   298.0,
}

for i in os.listdir("."):
    n_steps = 50000
    if i in TEMP.keys():
        filename = os.path.join(i, 'run.mdp')
        o = open(filename, 'w')
        lines = MDP.replace('%TEMP%', '%6.1f'%TEMP[i])
        lines = lines.replace('%STEPS%', '%d'%n_steps)
        o.write(lines)


