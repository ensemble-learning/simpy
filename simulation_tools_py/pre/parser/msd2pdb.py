from msd import parse_msd
from msd import sort_mol
from pdb import output_pdb

def msd2pdb(pbc, coords):
    LABLE = 'ATOM'
    for i in range(len(coords)):
        a = coords[i].split()
        atom_no = i + 1
        atom_name = a[1][0]
        res_no = int(a[9])
        res_name = a[10]
        x = float(a[6])
        y = float(a[7])
        z = float(a[8])
        coords[i] = "%-6s%5d  %-3s %3s  %4d  %8.3f%8.3f%8.3f\n"%(LABLE,atom_no, atom_name, res_name, res_no, x, y, z)
    pbc = pbc.strip().replace(r'PBC:', 'CRYST1') + ' P1\n' 
    return pbc
    
if __name__ == '__main__':
    pbc, coords, connect = parse_msd('ice_cub_633.msd')
    sort_mol(coords)
    pbc = msd2pdb(pbc, coords)
    output_pdb(pbc, coords)
    
    