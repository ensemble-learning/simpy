from ase.cluster.icosahedron import Icosahedron
from ase.io import write

ns = 36
atoms = Icosahedron('Cu', noshells=ns)

print('#atoms = {}'.format(len(atoms)))
write('icosahedron-%03d.xyz'%ns, atoms)

