from ase.cluster.decahedron import Decahedron
from ase.io import write

p = 15 # natoms on 100 face normal to 5-fold axis
q = 1  # natoms 0n 100 parallel to 5-fold axis
r = 0  # depth of the Marks re-entrance?
atoms = Decahedron('Cu', p, q, r)  
print('#atoms = {}'.format(len(atoms)))
#write('images/decahedron.png', atoms)
write('decahedron-%03d.xyz'%p, atoms)

"""
p = 15 6 nm  2815
p = 45 20 nm 75945
p = 80 32 nm 426680
"""
