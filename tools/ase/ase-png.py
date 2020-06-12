#@note: https://wiki.fysik.dtu.dk/ase/dev/_modules/ase/io/pov.html
#@note: https://wiki.fysik.dtu.dk/ase/ase/io/io.html
#@note: https://gitlab.com/ase/ase/-/merge_requests/1285

from ase import Atoms
from ase.io import read, write
from ase.io.pov import get_bondpairs, set_high_bondorder_pairs

struct = read('CONTCAR')
r = [{'O': 0.4, 'H':0.2, 'Cl':0.6}[at.symbol] for at in struct]

bondpairs = get_bondpairs(struct, radius=0.80)
#high_bondorder_pairs ={}
#high_bondorder_pairs[(0, 1)] = ((0, 0, 0), 2, (0.17, 0.17, 0))
#bondpairs = set_high_bondorder_pairs(bondpairs, high_bondorder_pairs)

#write('output.png', struct, **kwargs)

kwargs = {
    'rotation': '10z,-80x',
    'show_unit_cell':1,
    'radii':r,
    'colors':None,
    'bondatoms': bondpairs,
    'canvas_width': 1000, 
    }
kwargs.update({
    'run_povray': True,
    'display': False,
    'transparent': False,
    'background': 'White',
    'camera_dist': 10,
    })

write('o2.pov', struct*(1,1,1), **kwargs)
