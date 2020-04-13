import os
import numpy as np
from pgroup.ase.io import POV
from ase.build import fcc111

atoms = fcc111('Pt', (3, 2, 2), vacuum=10.)
atoms[11].symbol = 'Au'

# Repeat, rotate and translate to desired position.
atoms = atoms.repeat((6, 8, 1))
atoms.rotate('z', np.pi, rotate_cell=True)
atoms.translate(-atoms[190].position)

# Do the ray-trace.
pov = POV(atoms,
          cameralocation=(0., -15., 10.),
          area_light=[(2., -25., 20.), 'White', 1.7, 1.7, 3, 3],
          tex='vmd',
          clipplane='y, -0.00',
          pixelwidth=640,
          look_at=(0., 0., 2.),
          )

pov.write('cut-alloy.png')
