import os
import os.path
import shutil

for i in os.listdir("./"):
    if os.path.isdir(i):
        shutil.copyfile("%s/tatb.gin"%i, "%s/tatb_back.gin"%i)
        f = open("%s/tatb_back.gin"%i, 'r')
        o = open("%s/tatb.gin"%i, 'w')

        for j in f:
            if j.startswith("C1        O4            0.2020    3.5600"):
                o.write("C1        O4        0.260903     2.90043   0.0   8.0 0.0 0.0\n")
            else:
                o.write(j)
        o.close()
        f.close()
