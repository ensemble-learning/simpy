timestep	2.0
fix             shakeH all shake 0.0001 20 500 m shakeOpts
print                          .
print ================================================
print "NPT dynamics with an isotropic pressure of 1atm."
print ================================================
print                       .

fix             2 all npt temp 300.0 300.0 100.0 iso 1.0 1.0 2000.0
thermo          100
thermo_style    multi
dump            1 all custom 1000 ${sname}.npt.lammps id type xu yu zu vx vy vz
run		2500000 # run for 5 ns
restart         50000 ${sname}.*.restart
run             2500000 # run for 5 ns
unfix           2
undump          1
unfix           shakeH
