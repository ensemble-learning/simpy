from pdb import parse_pdb
from pdb import format_msdcoord
from pdb import output_pdb

def drag_pdb(coords, step):
    coords_drag = []
    for i in range(len(coords)):
        tokens = coords[i].split()
        tokens[7] = float(tokens[7]) + step
        coords_drag.append(tokens)
    return coords_drag

def pdb2gro(coords, pbc):
    import math
    
    coords_gro = []

    for i in coords:
        a = i.strip().split()
        atom_no = int(a[1])
        mol = a[3]
        atom = a[2]
        mol_no = int(a[4])
        atom_x = float(a[5])/10 #nm
        atom_y = float(a[6])/10
        atom_z = float(a[7])/10
        coords_gro.append([atom_no, mol, atom, mol_no, atom_x, atom_y, atom_z])
        
    [a, b, c, alpha, beta, gamma ]= [float(i) for i in pbc.split()[1:7]]
    a = a / 10
    b = b / 10
    c = c / 10
    cosA = math.cos(alpha*math.pi/180)
    cosB = math.cos(beta*math.pi/180)
    cosG = math.cos(gamma*math.pi/180)
    sinG = math.sin(gamma*math.pi/180)
    ax = a 
    bx = b * cosG 
    by = b * sinG 
    cx = c * cosB 
    cy = c * (cosA - cosB * cosG) / sinG
    cz = math.sqrt(c*c - cx*cx - cy*cy)
    
    box = [ax, by, cz, 0, 0, bx, 0, cx, cy ]
    print 'ax', ax
    print 'bx', bx
    print 'by', by
    print 'cx', cx
    print 'cy', cy
    print 'cz', cz
    
    return coords_gro, box
    
def pdb2data(coords, pbc):#not complete
    o = open('out.data', 'w')
    head = []
    head.append("# CT generated LAMMPS data file\n\n")
    head.append("%d atoms\n"%len(coords))
    head.append("0 bonds\n0 angles\n0 dihedrals\n0 impropers\n\n") #need improve
    head.append("2  atom types\n0  bond types\n0  angle types\n0  dihedral types\n0  improper types\n\n") #need improve
    o.close()

def main():
    pbc, coords = parse_pdb('layer.pdb')
    coords_gro, box = pdb2gro(coords, pbc)
    from gro import output_gro
    output_gro(coords_gro, box)

if __name__ == '__main__':
    main()
    
    