python ~/soft/simpy/tools/vasp/get_potcar.py
python ~/soft/simpy/tools/vasp/get_kpoint.py
python ~/soft/simpy/tools/vasp/incar.py

# for vasp 5.2
#genpbs_vasp

# for vasp 5.3
python /project/source/VASP/vasp.5.3.5/bin/genVASPpbs.py 1 4 24:00:00 vasp
