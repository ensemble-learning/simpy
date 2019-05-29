import shutil

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
    shutil.copy(flist[i], flist[i]+".0")
    o = open(flist[i], "w")
    for k in dtrain:
        o.write(k)
    for k in dtest:
        o.write(k)
    o.close()
    
        
    
    
    



