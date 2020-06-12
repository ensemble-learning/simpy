import os, sys

fs2ps = 0.001
n_atoms = 0
n_points = 0
# http://dospt.org/index.php/DoSPT

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
                    
# read input file

if not os.path.exists("inp"): 
    cells = get_cell()
    ts = get_time_step()
    tempt = get_temperature()
    o = open("inp", "w")
    o.write("cell =   %.3f %.3f %.3f # in nm\n"
            %(float(cells[0][0])/10.0, float(cells[1][1])/10.0, float(cells[2][2])/10.0))
    o.write("temperature = %.3f\n"%tempt)
    o.write("time_step = %.3f fs\n"%ts)
    o.write("selected-atoms = 196 49 98\n")
    sys.stderr.write("Fail to locate inp file, automatically generate a template.\n")
    sys.stderr.write("Please check the parameters\n")
                        
    sys.exit(0)
    
f = open("inp", "r")
lines = f.readlines()
f.close()
cell = " ".join(lines[0].strip().split()[2:5])
temperature = float(lines[1].strip().split()[2])
time_step = float(lines[2].strip().split()[2])
selected_atoms =  [int(i) for i in lines[3].strip().split()[2:]]

f = open("traj.xyz", "r")
n = 0
for i in f:
    if n == 0:
        n_atoms = int(i.strip())
        print(n_atoms)
    n += 1
n_points = int(n/(n_atoms+2))
f.close()

# write input file with comments
time_sim = n_points*fs2ps*time_step
o = open("input", "w")
o.write("# Number of points in the trajectory\n")
o.write("points = %d\n"%n_points)
o.write("# Trajectory period in ps\n")
o.write("tau = %.4f\n"%time_sim)
o.write("# Size of the box in nm\n")
o.write("cell = %s\n"%cell)
o.write("# Temperature in K\n")
o.write("temperature = %f\n"%temperature)
o.write("# Trajectory info\n")
o.write("format = xyz\n")
o.write("# Estimate velocities from positions\n")
o.write("estimate_velocities = .true.\n")

o = open("masses", "w")
o.write("""Cu  63.546
H   1.0080
Li  6.9400
C   12.011
O   15.999
F   18.998
Na  22.990
K   39.098
S   32.060
Au  196.97
N   14.007
Pt  195.08
Cl   35.45
Fe  55.845
B   10.81
""")
o.close()

o = open("supergroups", "w")
o.write("1\n")
o.close()

n_selected = len(selected_atoms)
o = open("groups", "w")
o.write("%d 2\n"%n_atoms)
o.write("%d 1\n"%n_selected)
for i in selected_atoms:
    o.write("%d "%i)
o.write("\n")

o.write("%d 1\n"%(n_atoms - n_selected))
for i in range(1, n_atoms+1):
    if i not in selected_atoms:
        o.write("%d "%i)
o.write("\n")
o.close()

o = open("run.sh", "w")
o.write("DoSPT\n")
o.write("python2 ~/soft/simpy/tools/dos2pt/plot_dos.py\n")
o.close()

os.system("chmod +x run.sh")
os.system("chmod +x ./run.sh")

