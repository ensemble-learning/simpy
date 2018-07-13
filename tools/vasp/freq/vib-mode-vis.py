import os
import numpy as np

os.system("xdatgen")

def get_atom_info():
    f = open("POSCAR", "r")
    lines = f.readlines()
    f.close()
    s = float(lines[1])
    a = np.array([float(i) for i in lines[2].strip().split()])
    b = np.array([float(i) for i in lines[3].strip().split()])
    c = np.array([float(i) for i in lines[4].strip().split()])
    a = a*s
    b = b*s
    c = c*s
    atom_lists = []
    ele = lines[5].strip().split()
    n_atoms = [int(i) for i in lines[6].strip().split()]
    for i in range(len(ele)):
        for j in range(n_atoms[i]*4):
            atom_lists.append(ele[i])
    return a, b, c, atom_lists

def read_frames(xdat_file):
    frames = []
    f = open(xdat_file, "r")
    lines = f.readlines()
    f.close()
    frame = []
    for i in lines[6:]:
        tokens = [float(j) for j in i.strip().split()]
        if len(tokens) == 0:
            frame = []
            frames.append(frame)
        else:
            frame.append(tokens)
    return frames

def write_modes(frames, a, b, c, atom_lists, mode_file):
    o = open(mode_file, "w")
    for i in range(len(frames)-1):
        o.write("%d\n\n"%len(atom_lists))
        for j in range(len(frames[i])):
            x = frames[i][j][0]
            y = frames[i][j][1]
            z = frames[i][j][2]
            x = a[0]*x + b[0]*y + c[0]*z
            y = a[1]*x + b[1]*y + c[1]*z
            z = a[2]*x + b[2]*y + c[2]*z
            o.write("%-6s"%atom_lists[j])
            o.write("%12.4f%12.4f%12.4f\n"%(x, y, z))
    o.close()

def main():
    # read pbc and atom lists
    a, b, c, atom_lists = get_atom_info()
    n_mode = 1

    while(1):
        xdat_file = "%03d.XDATCAR"%n_mode
        mode_file = "mode-%03d.xyz"%n_mode
        if os.path.exists(xdat_file):
            # read frames
            frames = read_frames(xdat_file)
            # output vib modes
            write_modes(frames, a, b, c, atom_lists, mode_file)
            n_mode += 1
        else:
            break

if __name__ == "__main__":
    main()

