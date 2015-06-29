"""
Generate DFTB restart files for new calculation.
"""
import os
import shutil

def get_last_vels():
    vels = []
    natom = 0
    f = open("geo_final.xyz", "r")

    n = 0
    for i in f:
        if n == 0:
            natom = int(i.strip())
        else:
            tokens = i.strip().split()
            if len(tokens) == 8:
                vx = float(tokens[5])
                vy = float(tokens[6])
                vz = float(tokens[7])
                vels.append([vx, vy, vz])
        n += 1

    f.close()

    o = open("vels", "w")
    for i in vels[-natom:]:
        o.write("%16.9f%16.9f%16.9f\n"%(i[0], i[1], i[2]))
    o.close()

def gen_restart():
    if not os.path.exists("restart"):
        os.mkdir("restart")
    path = os.path.join(os.getcwd(), "restart")
    shutil.copy("geo_end.gen", os.path.join(path, "geo_start.gen"))
    shutil.copy("charges.bin", os.path.join(path, "charges.bin"))
    get_last_vels()
    shutil.copy("vels", os.path.join(path, "vels"))

def main():
    gen_restart()

if __name__ == "__main__":
    main()
