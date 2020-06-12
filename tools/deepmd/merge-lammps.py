import shutil, os

folders = []
f = open("flist", "r")
for i in f:
    folders.append(i.strip())
f.close()

raw_files = ["energy.raw", "force.raw", "virial.raw"]
for i in raw_files:
    lines = []
    print(folders)
    for j in folders:
        print(j)
        f = open("./%s/lammps_raw/%s"%(j, i), "r")
        for k in f:
            lines.append(k)
        f.close()
    o = open(i, "w")
    for k in lines:
        print(k)
        o.write(k)
    o.close()


