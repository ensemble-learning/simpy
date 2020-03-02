import ase.io.vasp

struct = ase.io.vasp.read_vasp('002-mp-126-sur-111.vasp')
ase.io.vasp.write_vasp('POSCAR', struct*(2,2,1),direct=True, vasp5=True)

