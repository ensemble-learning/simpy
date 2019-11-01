# generate vasp raw
cmdall "cp CONTCAR POSCAR"
cmdall "python ~/soft/simpy/tools/deepmd/vasp2raw.py"

# generate lammps raw
cmdall "python ~/soft/simpy/lib/e_2_contcar.py POSCAR"
copyall ffield control.reaxc lammps_input 
cmdall '/central/home/tcheng/soft/lammps/lammps-16Mar18/src/lmp_intel_cpu_intelmpi -in lammps_input'
cmdall "python ~/soft/simpy/tools/deepmd/reaxff-ml/lammps2raw.py"

# generate lv raw
cmdall "python ~/soft/simpy/tools/deepmd/reaxff-ml/lvdiff_real.py "

# gether data
python ~/soft/simpy/tools/deepmd/merge-sp.py
