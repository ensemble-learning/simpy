rm ./distance ./rhombi_old ./rhombi ./surface_old ./surface

gfortran distances.f -o distance
gfortran rhombi_new.f -o rhombi_old
gfortran rhombi_pbc.f -o rhombi
gfortran surface.f -o surface_old
gfortran surface_pbc.f -o surface
