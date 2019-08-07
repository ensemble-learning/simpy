import os

# convert lammps unit real to metal
# this is for deepmdkit training consistency
kcal2ev = 4.3363e-2
b2ev = 6.2415e-9
n_dump_head = 9

def parse_block(lines):
    forces = []
    for i in range(n_dump_head, len(lines)):
        tokens = lines[i].strip().split()
        for j in tokens[-3:]:
            forces.append(float(j))
    return forces

def get_n_atoms():
    dump_file = "dump.lammpstrj"
    stop_sign = "ITEM: NUMBER OF ATOMS"

    n_atoms = 0
    f = open(dump_file, "r")
    for i in f:
        if i.strip().startswith(stop_sign):
            break
    for i in f:
        tokens = i.strip().split()
        n_atoms = int(tokens[0])
        break
    f.close()
    return n_atoms

def get_dump_forces(n_atoms):
    dump_file = "dump.lammpstrj"
    forces_frames = []
    f = open(dump_file, "r")
    flag = 1
    lines = []
    n = 0
    while(flag):
        flag = 0
        for i in f:
            flag = 1
            lines.append(i)
            n += 1
            if n % (n_dump_head + n_atoms) == 0:
                forces = parse_block(lines)
                forces_frames.append(forces)
                n = 0 
                lines = []
                break
    f.close()
    return forces_frames

def get_log_pe_press():
    log_file = "log.lammps"
    stop_sign_1 = "PotEng"
    stop_sign_2 = "Loop time"

    energy_frames = []
    stress_frames = []

    f = open(log_file, "r")
    for i in f:
        if i.strip().startswith(stop_sign_1):
            break
    for i in f:
        if i.strip().startswith(stop_sign_2):
            break
        tokens = i.strip().split()
        energy_frames.append(float(tokens[0]))
        t2 = []
        for ii in tokens[1:]:
            t2.append(float(ii))
        stress_frames.append(t2)
    f.close()

    return energy_frames, stress_frames

def output_raw(forces_frames, energy_frames, stress_frames):
    output_folder = "lammps_raw"
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    os.chdir(output_folder)

    o = open("energy.raw", "w")
    for i in energy_frames:
        tokens = i*kcal2ev
        o.write("%.10f "%tokens+"\n")
    o.close()

    o = open("force.raw", "w")
    for i in forces_frames:
        for j in i:
            tokens = j*kcal2ev
            o.write("%.10f "%tokens)
        o.write("\n")
    o.close()

    o = open("virial.raw", "w")
    for i in stress_frames:
        xx, yy, zz, xy, xz, yz = [ii*b2ev for ii in i]
        yx, zx, zy = xy, xz, yz
        o.write("%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n"%(xx, xy, xz, yx, yy, yz, zx, zy, zz))
    o.close()
    os.chdir("..")

def main():
    n_atoms = get_n_atoms()
    forces_frames = get_dump_forces(n_atoms)
    energy_frames, stress_frames = get_log_pe_press()
    output_raw(forces_frames, energy_frames, stress_frames)

if __name__ == "__main__":
    main()


