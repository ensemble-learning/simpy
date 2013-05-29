#    1SOL     OW    1   0.324   2.972   9.872

def parse_gro(file):
    """This function is used to abstract the atom information from gro
    file. It reads atom number, molecular, atom name, molecular number
    and coordinate following the format of gro file. It does not read
    velocities information"""

    f = open(file, 'r')
    lines = f.readlines()
    f.close()
    
    coords = []

    for i in range(2,len(lines)-1 ):
        atom_no = int(lines[i][:5])
        mol = lines[i][5:10].strip()
        atom = lines[i][10:15].strip()
        mol_no = int(lines[i][15:20])
        atom_x = float(lines[i][20:28])
        atom_y = float(lines[i][28:36])
        atom_z = float(lines[i][36:44])
        coords.append([atom_no, mol, atom, mol_no, atom_x, atom_y, atom_z])
     
    box = [float(i) for i in lines[-1].strip().split()]

    return coords, box
            
def output_gro(coords, box):
    o = open('out.gro', 'w')
    o.write('title\n') #may be insert time
    o.write('%5d\n'%len(coords))
    for i in range(len(coords)):
        o.write('%5d%5s%5s%5d%8.3f%8.3f%8.3f\n'%tuple(coords[i]))
    if len(box) == 3:
        o.write('%10.5f%10.5f%10.5f\n'%tuple(box))
    elif len(box) == 9:
        o.write('%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n'%tuple(box))
    o.close()

def add_wall(box, coords, r, pos, moltype, atomtype):
    add = []
    step_x = box[0] / int(box[0]/r)
    step_y = box[1] / int(box[1]/r)
    print step_x
    atom_no = coords[-1][0] + 1
    mol_no = coords[-1][3] + 1
    counter = 0

    for i in range(int(box[0]/r)):
        for j in range(int(box[1]/r)):
            x = step_x * i
            y = step_y * j
            z = pos 
            add_atom = [atom_no +counter , moltype, atomtype, mol_no, x, y,z]
            add.append(add_atom)
            counter +=1

    return add
    

if __name__ == '__main__':
    coords, box = parse_gro('em.gro')
    wall1 = add_wall(box, coords, 0.20, 3, 'W1', 'Na')
    coords = coords + wall1
    wall2 = add_wall(box, coords, 0.20, 12, 'W2', 'Li')
    coords = coords + wall2
    for i in range(1000):
        coords[i][-1] = coords[i][-1] - 2.1
    output_gro(coords, box)
