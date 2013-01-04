import sys
import os

counter = 0
normal_list = []

o = open("fit.log", 'w')
p = os.popen('tail -n 2 ./a*/result.log | grep -B 1 Job', 'r')
for i in p:
    if i.startswith("==>"):
        normal_list.append(i.split("/")[1])
p.close()

for i in normal_list:
    f = open("./%s/result.log"%i)
    o.write("%s\n"%i)
    for j in f:
        if j.strip().startswith("Final cell parameters and derivatives :"):
            counter = 2
        if counter > 0:
            if j.strip().startswith("----"):
                counter = counter - 1
            o.write(j)
    
    f.close()

o.close()
