import sys
import numpy as np
from ase.io import read

def draw_lines(outfile, atoms):
    data = np.loadtxt('lines.dat', dtype='int')
    for i in data:
        x0 = atoms[i[0]].position
        x1 = atoms[i[1]].position
        outfile.write('draw line {%.4f %.4f %.4f} {%.4f %.4f %.4f} width 3 style dashed\n'%(*x0, *x1))

def render_tachyon(outfile):
    outfile.write('render Tachyon scene.dat tachyon -aasamples 12 %s -format png -res 1800 1800 -o s18.png\n')
    outfile.write('render Tachyon scene.dat tachyon -aasamples 12 %s -format png -res 900 900 -o s09.png\n')

def main():
    fname = sys.argv[1]
    atoms = read(fname)

    with open('vmd-add.vmd', 'w') as outfile:
        draw_lines(outfile, atoms)
        render_tachyon(outfile)

if __name__ == '__main__':
    main()
