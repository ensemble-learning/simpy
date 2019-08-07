import os, shutil

def new_force(f0, f1):
    o = open("./lv_raw/force.raw", "w")
    nframe = 0
    if len(f0) == len(f1):
        nframe = len(f0)
    for i in range(nframe):
        t1 = f0[i].strip().split()
        t2 = f1[i].strip().split()
        nf = 0
        if len(t1) == len(t2):
            nf = len(t1)
        for j in range(nf):
            new_f = float(t1[j]) - float(t2[j])
            o.write("%.10f "%new_f)
        o.write("\n")
    o.close()

def new_energy(f0, f1):
    o = open("./lv_raw/energy.raw", "w")
    if len(f0) == len(f1):
        nframe = len(f0)
    for i in range(nframe):
        t1 = float(f0[i].strip())
        t2 = float(f1[i].strip())
        new_pe = t1 - t2
        o.write("%.10f\n"%new_pe)
    o.close()

def new_viral(f0, f1):
    o = open("./lv_raw/virial.raw", "w")
    nframe = 0
    if len(f0) == len(f1):
        nframe = len(f0)
    for i in range(nframe):
        t1 = f0[i].strip().split()
        t2 = f1[i].strip().split()
        nf = 0
        if len(t1) == len(t2):
            nf = len(t1)
        for j in range(nf):
            new_f = float(t1[j]) - float(t2[j])
            o.write("%.10f "%new_f)
        o.write("\n")
    o.close()

flag = 0

if os.path.exists("lammps_raw"):
    flag += 1

if os.path.exists("vasp_raw"):
    flag += 1

if flag == 2:
    if not os.path.exists("lv_raw"):
        os.mkdir("lv_raw")

    shutil.copy("./vasp_raw/box.raw", "./lv_raw")
    shutil.copy("./vasp_raw/coord.raw", "./lv_raw")
    shutil.copy("./vasp_raw/type.raw", "./lv_raw")

    f = open("./vasp_raw/force.raw", "r")
    force_vasp = f.readlines()
    f.close()

    f = open("./vasp_raw/energy.raw", "r")
    energy_vasp = f.readlines()
    f.close()

    f = open("./vasp_raw/virial.raw", "r")
    viral_vasp = f.readlines()
    f.close()

    f = open("./lammps_raw/force.raw", "r")
    force_lammps = f.readlines()
    f.close()

    f = open("./lammps_raw/energy.raw", "r")
    energy_lammps = f.readlines()
    f.close()

    f = open("./lammps_raw/virial.raw", "r")
    viral_lammps = f.readlines()
    f.close()

    new_force(force_vasp, force_lammps)
    new_energy(energy_vasp, energy_lammps)
    new_viral(viral_vasp, viral_lammps)

