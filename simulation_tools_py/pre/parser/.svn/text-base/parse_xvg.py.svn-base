def parse_xvg_traj(filename):
    #This is a code temperally handling coord readed from g_traj
    #format: parse_xvg_traj(filename) such as parse_xvg_traj("coord.xvg")
    coords = []
    f = open(filename, 'r')
    for i in f:
        if i.startswith("#") or i.startswith("@"):
            pass
        else:
            coords.append(i)
    return coords

def get_center(mol):
    sum = 0
    for i in range(len(mol)):
        sum += mol[i]
    return sum/len(mol)

def xvg_density(coords, box, nlayer, n_atoms, nskip):
    """compute the number density distribution along one axis according
    the information in coord.xvg file.
    format xvg_density(coords, box, nlayer, n_atoms)
    coords: coordination from xvg
    box: box length (nm)
    nlayer: number of layers
    n_atoms: number of atoms in one molecular
    """

    step = box / nlayer
    layer = [0]*(nlayer+1)
    nmol = 0
    if (len(coords[-1].split())-1) % (n_atoms*nskip) == 0:
        nmol = (len(coords[-1].split())-1)/n_atoms/nskip
    else:
        print "Warning the atom number supplied may be wrong"
    
    for i in range(len(coords)):
        tokens = coords[i].split()
        for j in range(nmol):
            mol = []
            for k in range(n_atoms):
                mol.append(float(tokens[n_atoms*nskip*j + k + 1]))
            center = get_center(mol)               
            layer[int(center/step)] += 1

    for i in range(len(layer)):
        layer[i] = float(layer[i]) / len(coords)
  
    return layer

def coord_average_coord(coords):
    """compute the average coordinate of every particles. Center the system every
    frame.
    g_traj -s run.tpr -f traj.trr -n index.ndx -ox coord -noy
    format average_coord(coords)"""

    n = (len(coords[0].split()) - 1)/2
   
    for i in range(len(coords)):
        min = 500
        tokens = coords[i]
        layer = xvg_density(coords, 15.0, 150, 2)
        for j in range(len(layer)):
            if min > layer[j]:
                min = layer[j]
                index = j
        center_shift = j * 0.1
        for j in range(n):
            x = float(tokens.split()[2*j + 1])
            if x < center_shift*0.1:
                x += 15.0
            z = float(tokens.split()[2*j + 2])
            #print x, z
        print "------------------------"
    #for i in range(len(coords)):
    tokens = coords[-1].split()
    #center = 0_z[i]))
    o.close()

def average_coord(coords):
    """compute the average coordinate of every particles. Center the system every
    frame.
    g_traj -s run.tpr -f traj.trr -n index.ndx -ox coord -noy
    format average_coord(coords)"""
    grid_size = 0.0375
    l_x = 15.
    l_z = 9.0
    grid = [[0.0]*(int(l_z/grid_size)+1) for i in range(int(l_x/grid_size)+1)]
    
    n = (len(coords[0].split()) - 1)/2
    
    for i in range(len(coords)):
        min = 500
        tokens = coords[i]
        layer = xvg_density([tokens], 15.0, 150, 1, 2)

        atom_x = [0]*n
        atom_z = [0]*n
        
        for j in range(len(layer)):
            if min > layer[j]:
                min = layer[j]
                index = j
        print "---------step%6d---------"%i
        center_shift = index * 0.1
        center = 0
        for j in range(n):
            x = float(tokens.split()[2*j + 1])
            if x < center_shift:
                x += 15.0
            z = float(tokens.split()[2*j + 2])
            atom_x[j] = x
            atom_z[j] = z
            center += x
        center = center / n
        
        for j in range(n):
            atom_x[j] = atom_x[j] - center + 7.5
            if atom_x[j] < 0:
                atom_x[j] += 15
            elif atom_x[j] > 15:
                atom_x[j] = atom_x[j] - 15
            index_a = int(atom_x[j]/grid_size)
            index_b = int((atom_z[j]-3.0)/grid_size)
            grid[index_a][index_b] += 1
    
    return grid
def parse_xvg():
    """ a general function which can read date from xvg and regulate them into structure"""
    tips = []
    detail = []
    terms = []
    rows = []
    f = open('energy.xvg', 'r')
    counter = 0
    for i in f:
        if i.startswith('#'):
            tips.append(i)
        elif i.startswith('@'):
            if i.startswith('@ s'):
                counter += 1
                terms.append(i.split()[-1])
            detail.append(i)
        else:
            break
    for i in range(counter+1):
        rows.append([])
    for i in f:
        tokens = i.split()
        if len(rows) != len(tokens):
            return 0
        else:
            for j in range(len(rows)):
                rows[j].append(tokens[j])
    f.close()
    return terms, rows
    
def average_terms(terms, rows, step):
    o = open('ave_terms.txt', 'w')
    o.write('%14s'%'step' + ''.join(['%14s'%str(m).strip(r'"') for m in terms]) + '\n')
    a = len(rows)
    b = len(rows[0])
    counter = 0
    values = [0] * a
    
    for i in range(b):
        if counter % step == 0 and counter > 0:
                o.write(''.join(['%14.2f'%(float(m)/step) for m in values]) + '\n')
                values = [0] * a
        for j in range(a):
            values[j] += float(rows[j][i])
        counter += 1
    o.write(''.join(['%14.2f'%(float(m)/step) for m in values]) + '\n')

    o.close()
        
def main():
    terms, rows = parse_xvg()
    average_terms(terms, rows, 100)
    
if __name__ == "__main__":
    main()
    