"""
Get the bond, angle and torsion tables from data.lammps (from tod's code)
@todo: merge to simpy/lib/data.py
@note: debug case glycine
"""

import sys
import os

DEBUG_PATH = "../../testcode/glycine/"

def data2qmgeo():
    mol = sys.argv[1]
    
    atp_map = {}
    atom_map = []
    atom_map.append(0)
    
    f = open("data.lammps", "r")
    
    for i in f:
        if i.strip().startswith("Masses"):
            break
    
    for i in f:
        if i.strip().startswith("Pair Coeffs"):
            break
        tokens = i.strip().split()
        if len(tokens) == 0:
            pass
        else:
            atn = int(tokens[0])
            atp = tokens[3]
            atp_map[atn] = atp
    
    for i in f:
        if i.strip().startswith("Atoms"):
            break
        
    for i in f:
        if i.strip().startswith("Bonds"):
            break
        tokens = i.strip().split()
        if len(tokens) == 0:
            pass
        else:
            atp = int(tokens[2])
            atom_map.append(atp_map[atp])
    
    bonds = []
    for i in f:
        if i.strip().startswith("Angles"):
            break
        tokens = i.strip().split()
        if len(tokens) == 0:
            pass
        else:
            a1 = int(tokens[2]) 
            a2 = int(tokens[3])
            note = "# %s - %s"%(atom_map[a1], atom_map[a2])
            bonds.append("%-20s%6.2f%8d%8d %s\n"%(mol, 1.0, a1, a2, note))
    
    angles = []
    for i in f:
        if i.strip().startswith("Dihedrals"):
            break
        tokens = i.strip().split()
        if len(tokens) == 0:
            pass
        else:
            a1 = int(tokens[2])
            a2 = int(tokens[3])
            a3 = int(tokens[4])
            note = "# %s - %s - %s"%(atom_map[a1], atom_map[a2], atom_map[a3])
            angles.append("%-20s%6.2f%8d%8d%8d %s\n"%(mol, 1.0, a1, a2, a3, note))
    
    torsions = []
    for i in f:
        if i.strip().startswith("Impropers"):
            break
        tokens = i.strip().split()
        if len(tokens) == 0:
            pass
        else:
            a1 = int(tokens[2])
            a2 = int(tokens[3])
            a3 = int(tokens[4])
            a4 = int(tokens[5])
            note = "# %s - %s - %s - %s"%(atom_map[a1], atom_map[a2], atom_map[a3], atom_map[a4])
            torsions.append("%-20s%6.2f%8d%8d%8d%8d %s\n"%(mol, 1.0, a1, a2, a3, a4, note))
    
    impro = []
    for i in f:
        tokens = i.strip().split()
        if len(tokens) == 0:
            pass
        else:
            a1 = int(tokens[2])
            a2 = int(tokens[3])
            a3 = int(tokens[4])
            a4 = int(tokens[5])
            note = "# %s - %s - %s - %s"%(atom_map[a1], atom_map[a2], atom_map[a3], atom_map[a4])
            impro.append("%-20s%6.2f%8d%8d%8d%8d %s\n"%(mol, 1.0, a1, a2, a3, a4, note))
    
    o = open("qm_geo", "w")
    
    for i in bonds:
        o.write(i)
    
    for i in angles:
        o.write(i)
    
    for i in torsions:
        o.write(i)
    
    for i in impro:
        o.write(i)
    
    o.close()

def main():
    data2qmgeo()

def debug():
    os.chdir(DEBUG_PATH)

if __name__ == "__main__":
    data2qmgeo()

