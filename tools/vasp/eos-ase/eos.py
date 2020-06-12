import numpy as np
import os
from ase import Atoms
from ase.constraints import FixAtoms, FixBondLengths
from ase.md import MDLogger, Langevin
import ase.units as units
from ase.io.trajectory import Trajectory
from ase.io import write,read
from ase.io.vasp import read_vasp, write_vasp
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS, FIRE
from ase.visualize import view
from ase.constraints import Hookean
from ase.neb import NEB
from ase.optimize import MDMin
from ase.optimize import BFGS
from ase.calculators.emt import EMT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

os.environ['VASP_COMMAND'] = 'mpirun -np 24 /opt/sourcecoude/vasp5.4.4/vasp.5.4.4/bin/vasp_std'
os.environ['VASP_PP_PATH'] = '/opt/software/vasp.5.4.4'

qmcalc = Vasp(
              system = 'chenggroup',
              prec = 'L',
              ncore = 4,
              xc = 'PBE',
              ivdw = 12,
              ismear = 0,
              sigma = 0.1,
              algo = 'Fast',
              nelm = 1000,
              lreal = 'Auto',
              encut = 400.0,
              ediff = 0.1e-04,
              isif = 2,
              isym = 0,
              lwave = False,
              lcharg = False,
              kpts = (9,9,9),
              gamma= True,
              ispin = 1,
              )

atoms = read_vasp('POSCAR')
atoms.set_calculator(qmcalc)
cell = atoms.get_cell()

traj = Trajectory('eos.traj', 'w')
for x in np.linspace(0.90, 1.10, 11):
    atoms.set_cell(cell * x, scale_atoms=True)
    atoms.get_potential_energy()
    folder = 'r%04d'%(int(x*1000))
    if not os.path.exists(folder):
        os.mkdir(folder)
    os.system("cp * %s"%folder)
    traj.write(atoms)
