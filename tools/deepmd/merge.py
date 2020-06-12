import shutil, os

def merge_dpdata():
    folders = []
    f = open("flist", "r")
    for i in f:
        if 0:
            if '.' in i:
                folders.append(i.split('.')[0]) 
            else:
                folders.append(i.strip())
        if 1:
            folders.append(i.strip())
    f.close()

    if os.path.exists("./%s/type.raw"%folders[0]):
        shutil.copy("./%s/type.raw"%folders[0], ".")
    if os.path.exists("./%s/type_map.raw"%folders[0]):
        shutil.copy("./%s/type_map.raw"%folders[0], ".")

    raw_files = ["box.raw", "coord.raw", "energy.raw", "force.raw", "virial.raw"]
    for i in raw_files:
        lines = []
        #print(folders)
        for j in folders:
            #print(j)
            if os.path.exists("./%s/%s"%(j, i)):
                f = open("./%s/%s"%(j, i), "r")
                for k in f:
                    lines.append(k)
                f.close()
        if len(lines) > 0:
            o = open(i, "w")
            for k in lines:
                o.write(k)
            o.close()
    if not os.path.exists('../01-raw'):
        os.mkdir('../01-raw')
    os.system('cp *.raw ../01-raw')

def merge_this():
    folders = []
    f = open("flist", "r")
    for i in f:
        if 0:
            if '.' in i:
                folders.append(i.split('.')[0]) 
            else:
                folders.append(i.strip())
        if 1:
            folders.append(i.strip())
    f.close()

    if os.path.exists("./%s/vasp_raw/type.raw"%folders[0]):
        print(folders)
        shutil.copy("./%s/vasp_raw/type.raw"%folders[0], ".")
    if os.path.exists("./%s/vasp_raw/type_map.raw"%folders[0]):
        shutil.copy("./%s/vasp_raw/type_map.raw"%folders[0], ".")

    raw_files = ["box.raw", "coord.raw", "energy.raw", "force.raw", "virial.raw"]
    for i in raw_files:
        lines = []
        #print(folders)
        for j in folders:
            #print(j)
            if os.path.exists("./%s/vasp_raw/%s"%(j, i)):
                f = open("./%s/vasp_raw/%s"%(j, i), "r")
                for k in f:
                    lines.append(k)
                f.close()
        if len(lines) > 0:
            o = open(i, "w")
            for k in lines:
                o.write(k)
            o.close()
    if not os.path.exists('../vasp_raw'):
        os.mkdir('../vasp_raw')
    os.system('cp *.raw ../vasp_raw')

merge_this()
