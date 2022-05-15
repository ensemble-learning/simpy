import os
from ase.io import read

fs2ps = 0.001

def get_cell():
    cells = []
    if os.path.exists('POSCAR'):
        f = open("POSCAR", "r")
        lines = f.readlines()
        f.close()
        for i in range(2,5):
            cells.append(lines[i].strip().split())
    return cells

def get_time_step():
    if os.path.exists('INCAR'):
        f = open("INCAR", "r")
        for i in f:
            if "POTIM" in i:
                tokens = i.strip().split("=")
                tokens = tokens[1].split()
                ts = float(tokens[0])
        f.close()
    else:
        ts = 5
    return ts

def get_temperature():
    if os.path.exists('INCAR'):
        f = open("INCAR", "r")
        for i in f:
            if "TEBEG" in i:
                tokens = i.strip().split("=")
                tokens = tokens[1].split()
                tempt = float(tokens[0])
        f.close()
    else:
        tempt = 300
    return tempt

def get_npoints(atoms):
    f = open("traj.xyz", "r")
    lines = f.readlines()
    f.close()
    n_points = int(len(lines)/(len(atoms)+2))
    return(n_points)

def gen_groups(atoms):
    os.system('python ~/soft/reactions/src/find-molecules.py POSCAR')
    mols = []
    if os.path.exists('mol.dat'):
        f = open('mol.dat', 'r')
        for i in f:
            tokens = i.strip().split()
            if len(tokens) > 2:
                mols.append(tokens)
        f.close()

    o = open('groups', 'w')
    o.write('%d %d\n'%(len(atoms), len(mols)))
    for n, i in enumerate(mols):
        mol_atoms = []
        for j in i[2:]:
            mol_atoms.append(int(j) + 1)
        o.write('%d 1\n'%len(mol_atoms))
        o.write(' '.join([str(jj) for jj in mol_atoms]) + '\n')
    o.close()

    o = open('supergroups', 'w')
    o.write(','.join([str(ii) for ii in range(1,len(mols)+1)])+'\n')
    o.close()

def gen_masses(atoms):
    elements, masses = [], []
    for i in atoms:
        if i.symbol not in elements:
            elements.append(i.symbol)
            masses.append(i.mass)

def gen_sh_run():
    o = open("run.sh", "w")
    o.write("DoSPT\n")
    o.write("/opt/software/anaconda3/bin/python ~/soft/simpy/tools/dos2pt/plot_dos.py\n")
    o.close()

    os.system("chmod +x run.sh")
    os.system("chmod +x ./run.sh")

def gen_input(atoms):
    # write input file with comments
    time_step = get_time_step()
    n_points = get_npoints(atoms)
    cell = atoms.get_cell()
    temperature = get_temperature()

    time_sim = n_points*fs2ps*time_step
    o = open("input", "w")
    o.write("# Number of points in the trajectory\n")
    o.write("points = %d\n"%n_points)
    o.write("# Trajectory period in ps\n")
    o.write("tau = %.4f\n"%time_sim)
    o.write("# Size of the box in nm\n")
    o.write("cell = %.2f %.2f %.2f\n"%(cell[0][0], cell[1][1], cell[2][2]))
    o.write("# Temperature in K\n")
    o.write("temperature = %f\n"%temperature)
    o.write("# Trajectory info\n")
    o.write("format = xyz\n")
    o.write("# Estimate velocities from positions\n")
    o.write("estimate_velocities = .true.\n")

atoms = read('POSCAR')
gen_groups(atoms)
gen_masses(atoms)
gen_input(atoms)
gen_sh_run()
