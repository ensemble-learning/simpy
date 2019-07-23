import shutil, os

flist = [ 
"box.raw",
"coord.raw",
"energy.raw",
"force.raw",
"virial.raw",
]

for i in range(len(flist)):
    f = open(flist[i], "r")
    data = f.readlines()
    f.close()

    dtrain, dtest = [], []

    for j in range(len(data)):
        if j%3 == 0:
            dtest.append(data[j])
        else:
            dtrain.append(data[j])
    # backup the files
    old_file = flist[i]
    n = 0
    flag = 1
    while(flag):
        new_file = old_file + ".%d"%n
        if not os.path.exists(new_file):
            flag = 0
        n += 1
    
    shutil.copy(old_file, new_file)
    o = open(flist[i], "w")
    for k in dtrain:
        o.write(k)
    for k in dtest:
        o.write(k)
    o.close()
    
        
    
    
    



