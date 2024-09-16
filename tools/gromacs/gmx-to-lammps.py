#!/usr/bin/env /opt/anaconda3/bin/python3

import os
from ase.io import read, write
from ase.data import atomic_masses, atomic_numbers

def write_in():
    o = open('in.sp', 'w')
    o.write('''units real
atom_style full

dimension 3
boundary p p p

bond_style hybrid harmonic
angle_style hybrid charmm
dihedral_style hybrid nharmonic
improper_style hybrid harmonic
special_bonds lj 0.0 0.0 0.5 coul 0.0 0.0 0.5

read_data lammps.data

pair_style lj/cut/coul/long 9.0 9.0
pair_modify tail yes
kspace_style pppm 1e-8

include      pair.in
include      bond-coeffs.in
include      angle-coeffs.in
include      dihedral-coeffs.in
include      improper-coeffs.in

pair_modify mix geometric

thermo_style custom ebond eangle edihed eimp epair evdwl ecoul elong etail pe

run 0
''')
    o.close()

def set_atps(atps):
    atps_sorted = []
    for i in atps:
        if not i in atps_sorted:
            atps.append(i)
    return atps_sorted

def n_to_mass(an):
    if an in atomic_numbers.values():
        return atomic_masses[an]
    else:
        return 'Invalid atomic number'

def parse_itp_file(file_path):
    # Initialize a dictionary to store the parsed data
    parsed_data = {}
    current_section = None

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            # Skip empty lines or comments
            if not line or line.startswith(';'):
                continue

            # Check for section headers
            if line.startswith('[') and line.endswith(']'):
                current_section = line[1:-1].strip()
                parsed_data[current_section] = []
            else:
                # Add data to the current section
                if current_section:
                    parsed_data[current_section].append(line)

    return parsed_data

atoms = read('conf.gro')

topol = parse_itp_file('topol.top')
atp = parse_itp_file('atp.itp')
atp_map = {}
for i in atp['atomtypes']:
    tokens = i.strip().split()
    if len(tokens) > 0:
        atp_map[tokens[0]] = tokens[1:]

mols = {}
for i in topol['defaults'][1:]:
    fname = i.strip().split()[1][1:-1]
    if fname != 'atp.itp':
        mol = parse_itp_file(fname)
        mol_name = mol['moleculetype'][0].split()[0]
        mols[mol_name] = mol

sim = []
for i in topol['molecules']:
    sim.append(i.strip().split())

mol_ids, atps, charges = [], [], []
cid = 1
for n, i in enumerate(sim):
    mol_name = i[0]
    n_mol = int(i[1])
    if n_mol > 0:
        for nn in range(n_mol):
            n_atoms = 0
            for ii in mols[mol_name]['atoms']:
                tokens = ii.strip().split()
                if len(tokens) > 0:
                    mol_ids.append(cid)
                    atps.append(tokens[1])
                    q = float(tokens[6])
                    charges.append(q)
                    n_atoms += 1
            cid += 1
        sim[n].append(n_atoms)
    else:
        # we just need to add a fick number to occupy the position
        sim[n].append(0)

atps_sorted = list(set(atps))
atps_sorted.sort()
atps_lammps = []
for i in atps:
    atps_lammps.append(atps_sorted.index(i)+1)

# write the pair coeffs
o = open('pair.in', 'w')
for n, i in enumerate(atps_sorted):
    o.write('pair_coeff    ')
    delta = float(atp_map[i][4])*10
    epsilon = float(atp_map[i][5])/4.184
    o.write('%4d %4d %12.8f %12.8f\n'%(n+1, n+1, epsilon, delta))
o.write('\n')
o.close()

# write the cell information
o = open('cell.data', 'w')
o.write('0.0 %12.7f xlo xhi\n'%atoms.cell[0][0])
o.write('0.0 %12.7f ylo yhi\n'%atoms.cell[1][1])
o.write('0.0 %12.7f zlo zhi\n'%atoms.cell[2][2])
if abs(atoms.cell[1][0]) > 0.01 or abs(atoms.cell[2][0]) > 0.01 or abs(atoms.cell[2][1]) > 0.01:
    o.write('%12.7f %12.7f %12.7f xy xz yz\n'%(atoms.cell[1][0], atoms.cell[2][0], atoms.cell[2][1]))
o.write('\n')
o.close()

# write the masses section
o = open('masses.data', 'w')
o.write('Masses\n\n')

for n, i in enumerate(atps_sorted):
    o.write('%9d '%(n+1))
    an = int(atp_map[i][0])
    mass = n_to_mass(an)
    o.write('%10.4f '%mass)
    o.write('   # %s'%i)
    o.write('\n')
o.write('\n')

o.close()
# write the atoms section
o = open('atoms.data', 'w')
o.write('Atoms # full\n\n')
for n, i in enumerate(atoms):
    o.write('%9d '%(n+1))
    o.write('%6d '%(mol_ids[n]))
    o.write('%4d '%(atps_lammps[n]))
    o.write('%7.3f '%(charges[n]))
    o.write('%24.14f '%i.x)
    o.write('%24.14f '%i.y)
    o.write('%24.14f '%i.z)
    o.write('\n')
o.write('\n')
o.close()

# get bond coeffs
bonds = []
for n, i in enumerate(sim):
    mol_name = i[0]
    n_mol = int(i[1])
    bonds.append([])
    if n_mol > 0:
        n_bonds = 0
        if 'bonds' in mols[mol_name].keys():
            for ii in mols[mol_name]['bonds']:
                tokens = ii.strip().split()
                if len(tokens) > 0:
                    bonds[n].append(tokens)
                    n_bonds += 1
        sim[n].append(n_bonds)
    else:
        # we just need to add a fake number to occupy the position
        sim[n].append(0)

# get the bonds
bond_type = 1
bond_types_lammps = []
for n, i in enumerate(bonds):
    if len(i) > 0:
        for ii in i:
            bl = float(ii[3])*10
            bk = float(ii[4])/4.184/2/100 # I dont know why 100 is needed here!
            bond_types_lammps.append([bond_type, 'harmonic', bk, bl])
            ii.append(bond_type)
            bond_type += 1

# write bond coeffs
o = open('bond-coeffs.in', 'w')
for i in bond_types_lammps:
    o.write('bond_coeff %4d  %s %12.4f %10.4f\n'%(i[0], i[1], i[2], i[3]))
o.write('\n')
o.close()

dn = 0
bonds_lammps = []
bond_id = 1
for n, i in enumerate(bonds):
    if len(i) > 0:
        for nn in range(int(sim[n][1])):
            for ii in i:
                bonds_lammps.append([bond_id, ii[-1], int(ii[0])+dn, int(ii[1])+dn])
                bond_id += 1
            dn += sim[n][2]
    else:
        dn += int(sim[n][1])*sim[n][2]

# write bonds
o = open('bonds.data', 'w')
if len(bonds_lammps) > 0:
    o.write('Bonds\n\n')
    for n, i in enumerate(bonds_lammps):
        o.write('%9d %9d %9d %9d\n'%(i[0], i[1], i[2], i[3]))
    o.write('\n')
o.close()

# get angle and angle coeffs
angles = []
for n, i in enumerate(sim):
    mol_name = i[0]
    n_mol = int(i[1])
    angles.append([])
    if n_mol > 0:
        n_angles = 0
        if 'angles' in mols[mol_name].keys():
            for ii in mols[mol_name]['angles']:
                tokens = ii.strip().split()
                if len(tokens) > 0:
                    angles[n].append(tokens)
                    n_angles += 1
        sim[n].append(n_angles)
    else:
        # we just need to add a fake number to occupy the position
        sim[n].append(0)

# get the angles 
angle_type = 1
angle_types_lammps = []
for n, i in enumerate(angles):
    if len(i) > 0:
        for ii in i:
            #['17', '7', '18', '5', '108.371', '256.916', '0.17891', '11828.774']
            theta = float(ii[4])
            k1 = float(ii[5])/4.184/2 
            r0 = float(ii[6])*10 
            k2 = float(ii[7])/4.184/2/100 # I dont know why 100 is needed here!
            angle_types_lammps.append([angle_type, 'charmm', k1, theta, k2, r0])
            ii.append(angle_type)
            angle_type += 1

# write angle coeffs
o = open('angle-coeffs.in', 'w')
for i in angle_types_lammps:
    o.write('angle_coeff %4d  %s %12.4f %10.4f %12.4f %10.4f\n'%(i[0], i[1], i[2], i[3], i[4], i[5]))
o.write('\n')
o.close()

dn = 0
angles_lammps = []
angle_id = 1
for n, i in enumerate(angles):
    if len(i) > 0:
        for nn in range(int(sim[n][1])): # 1 is the number of species
            for ii in i:
                angles_lammps.append([angle_id, ii[-1], int(ii[0])+dn, int(ii[1])+dn, int(ii[2])+dn])
                angle_id += 1
            dn += sim[n][2] # 2 is the number of atom of each specie
    else:
        dn += int(sim[n][1])*sim[n][2]

# write angles 
o = open('angles.data', 'w')
if len(angles_lammps) > 0:
    o.write('Angles\n\n')
    for n, i in enumerate(angles_lammps):
        o.write('%9d %9d %9d %9d %9d\n'%(i[0], i[1], i[2], i[3], i[4]))
    o.write('\n')
o.close()

# get dihedral and dihedral coeffs
dihedrals = []
for n, i in enumerate(sim):
    mol_name = i[0]
    n_mol = int(i[1])
    dihedrals.append([])
    if n_mol > 0:
        n_dihedrals = 0
        if 'dihedrals' in mols[mol_name].keys():
            for ii in mols[mol_name]['dihedrals']:
                tokens = ii.strip().split()
                if len(tokens) > 0:
                    if int(tokens[4]) == 3:
                        dihedrals[n].append(tokens)
                        n_dihedrals += 1
        sim[n].append(n_dihedrals)
    else:
        # we just need to add a fake number to occupy the position
        sim[n].append(0)

# get the dihedrals
dihedral_type = 1
dihedral_types_lammps = []
for n, i in enumerate(dihedrals):
    if len(i) > 0:
        for ii in i:
            #['17', '7', '18', '5', '108.371', '256.916', '0.17891', '11828.774']
            k1 = float(ii[5])/4.184 
            k2 = float(ii[6])/-4.184
            k3 = float(ii[7])/4.184
            k4 = float(ii[8])/-4.184
            k5 = float(ii[9])/4.184
            k6 = float(ii[10])/-4.184
            dihedral_types_lammps.append([dihedral_type, 'nharmonic', 6, k1, k2, k3, k4, k5, k6])
            ii.append(dihedral_type)
            dihedral_type += 1

# write dihedral coeffs
o = open('dihedral-coeffs.in', 'w')
for i in dihedral_types_lammps:
    #o.write('dihedral_coeff %4d  %s %12.8f %12.8f %12.8f %12.8f %12.8f\n'%(i[0], i[1], i[2], i[3], i[4], i[5], i[6]))
    o.write('dihedral_coeff %4d  %s %2d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n'%(i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8]))
o.write('\n')
o.close()

dn = 0
dihedrals_lammps = []
dihedral_id = 1
for n, i in enumerate(dihedrals):
    if len(i) > 0:
        for nn in range(int(sim[n][1])): # 1 is the number of species
            for ii in i:
                dihedrals_lammps.append([dihedral_id, ii[-1], int(ii[0])+dn, int(ii[1])+dn, int(ii[2])+dn, int(ii[3])+dn])
                dihedral_id += 1
            dn += sim[n][2] # 2 is the number of atom of each specie
    else:
        dn += int(sim[n][1])*sim[n][2]

# write dihedrals 
o = open('dihedrals.data', 'w')
if len(dihedrals_lammps) > 0:
    o.write('Dihedrals\n\n')
    for n, i in enumerate(dihedrals_lammps):
        o.write('%9d %9d %9d %9d %9d %9d\n'%(i[0], i[1], i[2], i[3], i[4], i[5]))
    o.write('\n')
o.close()

# get improper and improper coeffs
impropers = []
for n, i in enumerate(sim):
    mol_name = i[0]
    n_mol = int(i[1])
    impropers.append([])
    if n_mol > 0:
        n_impropers = 0
        if 'dihedrals' in mols[mol_name].keys():
            for ii in mols[mol_name]['dihedrals']:
                tokens = ii.strip().split()
                if len(tokens) > 0:
                    if int(tokens[4]) == 2:
                        impropers[n].append(tokens)
                        n_impropers += 1
        sim[n].append(n_impropers)
    else:
        # we just need to add a fake number to occupy the position
        sim[n].append(0)

# get the impropers
improper_type = 1
improper_types_lammps = []
for n, i in enumerate(impropers):
    if len(i) > 0:
        for ii in i:
            k1 = float(ii[5]) 
            k2 = float(ii[6])/4.184/2
            improper_types_lammps.append([improper_type, 'harmonic', k2, k1])
            ii.append(improper_type)
            improper_type += 1

# write improper coeffs
o = open('improper-coeffs.in', 'w')
for i in improper_types_lammps:
    o.write('improper_coeff %4d  %s %14.8f %10.4f\n'%(i[0], i[1], i[2], i[3]))
o.write('\n')
o.close()

dn = 0
impropers_lammps = []
improper_id = 1
for n, i in enumerate(impropers):
    if len(i) > 0:
        for nn in range(int(sim[n][1])): # 1 is the number of species
            for ii in i:
                impropers_lammps.append([improper_id, ii[-1], int(ii[0])+dn, int(ii[1])+dn, int(ii[2])+dn, int(ii[3])+dn])
                improper_id += 1
            dn += sim[n][2] # 2 is the number of atom of each specie
    else:
        dn += int(sim[n][1])*sim[n][2]

# write impropers 
o = open('impropers.data', 'w')
if len(impropers_lammps) > 0:
    o.write('Impropers\n\n')
    for n, i in enumerate(impropers_lammps):
        o.write('%9d %9d %9d %9d %9d %9d\n'%(i[0], i[1], i[2], i[3], i[4], i[5]))
    o.write('\n')
o.close()

# write the head section
o = open('head.data', 'w')
o.write('Au cas ou je ne vous verrais pas, bon apres-midi, bonne soiree, et bonne nuit.\n\n')
o.write('%d atoms\n'%len(atoms))
o.write('%d bonds\n'%len(bonds_lammps))
o.write('%d angles\n'%len(angles_lammps))
o.write('%d dihedrals\n'%len(dihedrals_lammps))
o.write('%d impropers\n'%len(impropers_lammps))
o.write('\n')
o.write('%d atom types\n'%len(atps_sorted))
o.write('%d bond types\n'%len(bond_types_lammps))
o.write('%d angle types\n'%len(angle_types_lammps))
o.write('%d dihedral types\n'%len(dihedral_types_lammps))
o.write('%d improper types\n'%len(improper_types_lammps))
o.write('\n')
o.close()

os.system('cat head.data cell.data masses.data atoms.data bonds.data angles.data dihedrals.data impropers.data > lammps.data')
write_in()
