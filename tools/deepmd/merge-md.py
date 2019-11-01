import shutil, os

raw_files = ["box.raw", "coord.raw", "energy.raw", "force.raw", "virial.raw"]

folders = []
f = open("flist", "r")
for i in f:
    folders.append(i.strip())
f.close()

for i in folders:
    os.chdir(i)
    os.chdir('poscars')
    data = []
    f = open("flist", "r")
    for j in f:
        data.append(j.strip())
    f.close()
    if not os.path.exists('vasp_raw'):
        os.mkdir('vasp_raw')
    if os.path.exists('./%sext/vasp_raw/type.raw'%data[0]):
        shutil.copy("./%sext/vasp_raw/type.raw"%data[0], "vasp_raw")
    for k in raw_files:
        lines = []
        for m in data:
            if os.path.exists("%sext/vasp_raw/%s"%(m, k)):
                f2 = open("%sext/vasp_raw/%s"%(m, k))
                for n in f2:
                    lines.append(n)
                f2.close()
        if len(lines) > 0:
            o = open('./vasp_raw/%s'%k, "w")
            for k in lines:
                o.write(k)
            o.close()

    if not os.path.exists('lv_raw'):
        os.mkdir('lv_raw')
    if os.path.exists('./%sext/vasp_raw/type.raw'%data[0]):
        shutil.copy("./%sext/vasp_raw/type.raw"%data[0], "lv_raw")

    for k in raw_files:
        lines = []
        for m in data:
            if os.path.exists("%sext/lv_raw/%s"%(m, k)):
                f2 = open("%sext/lv_raw/%s"%(m, k))
                for n in f2:
                    lines.append(n)
                f2.close()
        if len(lines) > 0:
            o = open('./lv_raw/%s'%k, "w")
            for k in lines:
                o.write(k)
            o.close()
    os.chdir('..')
    os.chdir('..')


if not os.path.exists('vasp_raw'):
    os.mkdir('vasp_raw')
if not os.path.exists('lv_raw'):
    os.mkdir('lv_raw')

if os.path.exists('./%s/poscars/vasp_raw/type.raw'%folders[0]):
    shutil.copy('./%s/poscars/vasp_raw/type.raw'%folders[0], 'vasp_raw')
    shutil.copy('./%s/poscars/vasp_raw/type.raw'%folders[0], 'lv_raw')

for i in raw_files:
    lines = []
    for j in folders: 
        if os.path.exists('%s/poscars/vasp_raw/%s'%(j, i)):
            f = open('%s/poscars/vasp_raw/%s'%(j, i), 'r')
            for k in f:
                lines.append(k)
            f.close()
    if len(lines) > 0:
        o = open('./vasp_raw/%s'%i, 'w')
        for k in lines:
            o.write(k)
        o.close()
            
    lines = []
    for j in folders: 
        if os.path.exists('%s/poscars/lv_raw/%s'%(j, i)):
            f = open('%s/poscars/lv_raw/%s'%(j, i), 'r')
            for k in f:
                lines.append(k)
            f.close()
    if len(lines) > 0:
        o = open('./lv_raw/%s'%i, 'w')
        for k in lines:
            o.write(k)
        o.close()
