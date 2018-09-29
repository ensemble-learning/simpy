def parse_pdb(file):
    """parse the msd file into three part: pbc, coordinates and 
    connectivity"""
    pbc = ''
    coords = []
    f = open(file, 'r')
    for i in f:
        if i.startswith('REMARK'):
            pass
        elif i.strip().startswith('CRYST'):
            pbc = i
        elif len(i.split()) == 11:
            coords.append(i)
        elif len(i.split()) == 10:
            coords.append(i)
        elif len(i.split()) == 8:
            coords.append(i)
    
    return pbc, coords
    
def output_pdb(pbc, coords):
    """output the msd file using pbc, coords and connect"""
    
    o = open('out.pdb', 'w')
    o.write(pbc)
    for i in coords:
        o.write(i)
    o.write('TER\n')
    o.close()

def format_msdcoord(coords):
    for i in range(len(coords)):
        #print coords
        lable = coords[i][0]
        atom_no = int(coords[i][1])
        atom_name = coords[i][2]
        res_name = coords[i][3]
        res_no = int(coords[i][4])
        x = float(coords[i][5])
        y = float(coords[i][6])
        z = float(coords[i][7])
        
        coords[i] ="%-6s%5d  %-3s %3s  %4d  %8.3f%8.3f%8.3f\n"%(lable,atom_no, atom_name, res_name, res_no, x, y, z)
            
if __name__ == '__main__':
    pbc, coords = parse_pdb('out.pdb')
    print pbc
    from spc_tip4p import spc_tip4p
    coords_tip4p = spc_tip4p(coords)
    format_msdcoord(coords_tip4p)
    output_pdb(pbc, coords_tip4p)
    
    
