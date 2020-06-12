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

os.environ['VASP_COMMAND'] = 'srun -n 128 /central/home/tcheng/soft/vasp.5.4.4/vasp.5.4.4/bin/vasp_gam'
os.environ['VASP_PP_PATH'] = '/central/home/tcheng/soft/vasp.5.3.5/vasp.5.3.5/potcar'

qmcalc = Vasp(
              system = 'chenggroup',
              prec = 'L',
              ncore = 4,
              xc = 'PBE',
              ivdw = 12,
              ismear = 0, # ismear = 1
              sigma = 0.1,
              algo = 'Fast',
              nelm = 1000,
              lreal = 'Auto',
              encut = 400.0, # encut
              ediff = 0.1e-04, # encut
              isif = 2,
              isym = 0,
              lwave = False,
              lcharg = False,
              kpts = (1,1,1), # kpoint
              gamma= True,
              ispin = 1,
              )

shutil.copy("POSCAR_0", "POSCAR")

tag = 'md'

# first run 1 ps NVT equilibration at 300 K
atoms = read_vasp('POSCAR')
atoms.set_calculator(qmcalc)

nvt = 1
if nvt:
    if not os.path.exists("nvt"):
        os.mkdir("nvt")
    os.chdir("nvt")
    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
    #md = nvtberendsen.NVTBerendsen(atoms, 1 * units.fs, 300, taut=0.5*1000*units.fs)
    md = Langevin(atoms, 1 * units.fs, temperature=300 * units.kB, friction=0.10)
    md.attach(MDLogger(md, atoms, tag+'.log', stress=True, mode='a'), interval=1)
    traj = Trajectory(tag + '.traj', 'w', atoms)
    md.attach(traj.write, interval=1)
    md.run(800)  # number of nvt steps 800
    shutil.copy("CONTCAR", "POSCAR")
    os.chdir("..")

atoms = read_vasp('POSCAR')
atoms.set_calculator(qmcalc)

o_flist = open("flist", "w")

for i in range(300, 900, 100):
    folder = "r%04d"%i
    if not os.path.exists(folder):
        os.mkdir(folder)
    o_flist.write("%s\n"%folder)
    os.chdir(folder)
    MaxwellBoltzmannDistribution(atoms, i * units.kB)
    #md = NVTBerendsen(atoms, 1 * units.fs, 300, taut=0.5*1000*units.fs)
    md = Langevin(atoms, 1 * units.fs, temperature=i * units.kB, friction=0.02)
    md.attach(MDLogger(md, atoms, tag+'.log', stress=True, mode='a'), interval=1)
    traj = Trajectory(tag + '.traj', 'w', atoms)
    md.attach(traj.write, interval=1)
    md.run(1000) # number of nvt steps at each temperature 1000

    atoms = read_vasp('CONTCAR')
    atoms.set_calculator(qmcalc)
    os.chdir("..")

o_flist.close()


