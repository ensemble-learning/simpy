import os,sys
from ase.io.trajectory import PickleTrajectory
from ase import Atoms
from ase.io.trajectory import Trajectory
import ase.io


#reading in the trajectory file created during dynamics
traj = Trajectory('%s'%sys.argv[1])

o = open("flist", "w")

ns = 0
for atoms in traj:
    fname = "ase_%04d.pdb"%ns
    o.write("%s\n"%fname)
    ase.io.write(fname, atoms)
    
    ns += 1
o.close()
