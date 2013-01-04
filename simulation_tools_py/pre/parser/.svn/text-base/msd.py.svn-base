def parse_msd(file):
    """parse the msd file into three part: pbc, coordinates and 
    connectivity"""
    pbc = ''
    coords = []
    connect = []
    f = open(file, 'r')
    for i in f:
        if i.startswith('#'):
            pass
        elif i.startswith('PBC'):
            pbc = i
        elif len(i.split()) == 12:
            coords.append(i)
        elif len(i.split()) == 3:
            connect.append(i)
    
    return pbc, coords, connect
def output_msd(pbc, coords, connect):
    """output the msd file using pbc, coords and connect"""
    
    o = open('out.msd', 'w')
    o.write("""#DFF:MSD
#Model Structure Data File    Energy = 0.0
""")
    o.write(pbc)
    o.write('%d\n'%len(coords))
    for i in coords:
        o.write(i)
    o.write('%d\n'%len(connect))
    for i in connect:
        o.write(i)
    o.write('M  END\n')
    o.close()
    
def sort_mol(coords):
    """ a method can sort the atom in msd according to the molecular
    it belong to and its atom type"""
    atom_quene = { 'O': 1, 'H': 2}
    for i in range(len(coords)):
        type = atom_quene[coords[i].split()[1]]
        coords[i] = ('%04d%s'%(int(coords[i].split()[-3]), type), coords[i])
    coords.sort()
    for i in range(len(coords)):
        coords[i] = coords[i][1]
        
if __name__ == '__main__':
    pbc, coords, connect = parse_msd('test.msd')
    sort_mol(coords)
    output_msd(pbc, coords, connect)
    
    