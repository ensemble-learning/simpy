import os

def auto_fit(atomType, parameter, scale_factor):
    itpFile = ""
    for i in scale_factor:
        if os.path.exists("%s_%03d"%(parameter, i)):
            pass
        else:
            os.mkdir("%s_%03d"%(parameter, i))
        os.system("cp run.mdp run.gro  *.top ccl4.itp na.itp %s_%03d"%(parameter, i))
    fileList = os.listdir("./")
    for i in fileList:
        if i[-3:] == "itp" and i.startswith("s"):
            itpFile = i
    if os.path.exists(itpFile):
        for i in scale_factor:
            f = open(itpFile, 'r')
            o = open("./%s_%03d/%s"%(parameter, i, itpFile), 'w')
            for j in f:
                if j.strip().startswith(atomType):
                    a = j.split()
                    if parameter == "delta":
                        nobond = float(a[4]) * (i / 100.0)
                        o.write("%8s%12s%12s%8s%12.4f%12s\n"%(a[0], a[1], a[2], a[3], nobond, a[5]))
                    elif parameter == "sigma":
                        nobond = float(a[5]) * (i / 100.0)
                        o.write("%8s%12s%12s%8s%12s%12.4f\n"%(a[0], a[1], a[2], a[3],a[4], nobond))
                elif j.strip().startswith("[ bondtypes ]") or j.strip().startswith("[ pairtypes ]"):
                    o.write(j)
                    break
                else:
                    o.write(j)
            for j in f:
                    o.write(j)
            o.close()
            f.close()

if __name__ == "__main__":
    auto_fit("c_2", "delta", [90, 110])
    auto_fit("c_2", "sigma", [70, 130])
