import sys
from ase.io import read
from scipy.constants import N_A

if len(sys.argv) == 4:
    atoms = read(sys.argv[1])
    rho = float(sys.argv[2])
    l = float(sys.argv[3])

    wt = 0.0
    for i in atoms.get_masses():
        wt += i

    l = l*1e-10
    v = l*l*l*1e6

    nmol = int(round(v*rho*N_A/wt))
    print(nmol)
else:
    print('Usage: estimate_nmol.py pdbfile rho (in g/cm3) l (in A)')

