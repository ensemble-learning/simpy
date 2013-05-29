import math
from gro import parse_gro
from gro import output_gro

def get_dis(p1, p2):
    x = p1[0] - p2[0]
    y = p1[1] - p2[1]
    z = p1[2] - p2[2]
    return math.sqrt(x*x + y*y + z*z)
    
def get_cent(mol):
    center = [0, 0, 0]
    for i in range(len(mol)):
        center[0] += mol[i][-3]
        center[1] += mol[i][-2]
        center[2] += mol[i][-1]
    center[0] = center[0]/len(mol)
    center[1] = center[1]/len(mol)
    center[2] = center[2]/len(mol)
    return center
    
    
def add_atoms(coords, mol, n, tor, boxh, boxl):
    """Divide the box into n*n*n grid according both to the number of added molecular.
    If there is other particles within the tolerance, the insertion fails."""
    
    #divide the space into grid
    grid = []
    grid_add = []
    atom_no = coords[-1][0]
    mol_no = coords[-1][0]
    n = int(math.pow(n, 1/3.0)) + 1
    x_in = float(boxh[0] -  boxl[0])/n
    y_in = float(boxh[1] -  boxl[1])/n
    z_in = float(boxh[2] -  boxl[2])/n
    if x_in < tor or y_in < tor or z_in < tor:
        print 'Error'
        print x_in, y_in, z_in
        return 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                x = i*x_in + boxl[0]
                y = j*y_in + boxl[1]
                z = k*z_in + boxl[2]
                grid.append([x, y, z])
    for i in range(len(grid)):
        counter = 0
        p1 = grid[i]
        for j in coords:
            p2 = j[-3:]
            if get_dis(p1, p2) < tor:
                counter += 1
        if counter == 0:
            grid_add.append(grid[i])
    # add atom
    center = get_cent(mol)
    for i in range(len(mol)):
        mol[i][-3] = mol[i][-3] - center[-3]
        mol[i][-2] = mol[i][-2] - center[-2]
        mol[i][-1] = mol[i][-1] - center[-1]

    for i in range(len(grid_add)):
        mol_no += 1
        for j in range(len(mol)):
            add = [1, 'ETH', 'C', 1, 0.0, 0.0, 0.0]
            atom_no += 1
            add[0] = atom_no
            add[1] = mol[j][1]
            add[2] = mol[j][2]
            add[3] = mol_no
            add[4] = mol[j][4] + grid_add[i][-3]
            add[5] = mol[j][5] + grid_add[i][-2]
            add[6] = mol[j][6] + grid_add[i][-1]
            
            coords.append(add)

def enlarge(coords, box, dx, dy, dz):
    for i in range(len(coords)):
        #coords.append([atom_no, mol, atom, mol_no, atom_x, atom_y, atom_z])
        coords[i][-3] = coords[i][-3] * dx
        coords[i][-2] = coords[i][-2] * dy
        coords[i][-1] = coords[i][-1] * dz
    box[0] = box[0] * dx
    box[1] = box[1] * dy
    box[2] = box[2] * dz
    
def main():
    coords, box = parse_gro('conf999.gro')
    enlarge(coords, box, dx, dy, dz):
    output_gro(coords, box)

if __name__ == '__main__':
    main()
    
    
