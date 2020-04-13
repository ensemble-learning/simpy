import os
import numpy as np
from ase import Atoms
from ase.io import read
from ase.optimize import BFGS, FIRE
from ase.vibrations import Vibrations
from ase.calculators.vasp import Vasp2
from ase.build import molecule

#from deepmd.calculator import DP

#os.environ['ASE_VASP_COMMAND'] = 'srun -n 4 /central/groups/wag/programs/vasp.5.4.4/bin/vasp_gam'
os.environ['ASE_VASP_COMMAND'] = '/central/groups/wag/programs/vasp.5.4.4/bin/vasp_std'
os.environ['VASP_PP_PATH'] = '/central/home/tcheng/soft/vasp.5.3.5/vasp.5.3.5/potcar'
vasp_folder = 'vasp-run'

def set_vasp_calc():
    # set up vasp calculations
    vasp_calc = Vasp2(kpts=(1, 1, 1), directory=vasp_folder, pp='pbe')
    vasp_calc.set(
                isym=0,           # Turn off kpoint symmetry reduction
                encut=400,
                npar=2,
                nelmin=5,
                ediff=1e-5,
                prec='N',
                nelm=60,
                ismear=0,
                sigma=0.1,
                ivdw=12,
                ispin=1,
                lreal='Auto',
                isif=2,
                gga='PE',
    )
    return vasp_calc

atoms = read('POSCAR')
calc = set_vasp_calc()

atoms.set_calculator(calc)
e0 = atoms.get_potential_energy()
print('e0: %.6f'%e0)
#dyn = FIRE(atoms)
dyn = BFGS(atoms)
dyn.run(fmax=0.001)
e1 = atoms.get_potential_energy()
print('e1: %.6f'%e1)
vib = Vibrations(atoms)
vib.run()
freqs = vib.get_frequencies()
o = open('freq.dat', 'w')
for i in range(6,len(freqs)):
    o.write('%.4f\n'%np.real(freqs[i]))
o.close()
vib.summary()
