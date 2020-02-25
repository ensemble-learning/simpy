import os
folder_name = os.path.basename(os.getcwd())
pe = []
f = open("OSZICAR", 'r')
for i in f:
    if "E0" in i:
        tokens = i.strip().split()
        for j in tokens:
            if 'E0' in j:
                e0 = tokens[tokens.index(j)+1]
                pe.append(float(e0))
f.close()

o = open('e0.dat', 'w')
o.write('%s %.6f\n'%(folder_name, pe[-1]))
o.close()
