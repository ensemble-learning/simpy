import shutil, os

raw_files = ["box.raw", "coord.raw", "energy.raw", "force.raw", "virial.raw"]

folders = []
f = open("flist", "r")
for i in f:
    folders.append(i.strip())
f.close()

if not os.path.exists('vasp_raw'):
    os.mkdir('vasp_raw')
if not os.path.exists('lv_raw'):
    os.mkdir('lv_raw')

if os.path.exists('./%s/vasp_raw/type.raw'%folders[0]):
    shutil.copy('./%s/vasp_raw/type.raw'%folders[0], 'vasp_raw')
    shutil.copy('./%s/vasp_raw/type.raw'%folders[0], 'lv_raw')

for i in raw_files:
    lines = []
    for j in folders: 
        if os.path.exists('%s/vasp_raw/%s'%(j, i)):
            f = open('%s/vasp_raw/%s'%(j, i), 'r')
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
        if os.path.exists('%s/lv_raw/%s'%(j, i)):
            f = open('%s/lv_raw/%s'%(j, i), 'r')
            for k in f:
                lines.append(k)
            f.close()
    if len(lines) > 0:
        o = open('./lv_raw/%s'%i, 'w')
        for k in lines:
            o.write(k)
        o.close()
