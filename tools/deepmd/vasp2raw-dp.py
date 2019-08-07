import dpdata

d_poscar = dpdata.System('POSCAR')
d_outcar = dpdata.LabeledSystem('OUTCAR')
d_outcar.to_deepmd_raw('vasp_raw')
