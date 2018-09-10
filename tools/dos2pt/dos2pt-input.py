import os, sys

fs2ps = 0.001
n_atoms = 0
n_points = 0
# http://dospt.org/index.php/DoSPT

# read input file

if not os.path.exists("inp"): 
    o = open("inp", "w")
    o.write("""cell =   1.02 1.02 4.00 # in nm
temperature = 298.00
time_step = 1.0 fs
selected-atoms = 196 49 98
""")
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
C   12.011
O   16.00
H   1.01
Na  23.00
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
