import os, math

n = 0
f = open("force.raw", 'r')
for i in f:
    if len(i.strip())> 0:
        n += 1
f.close()

n_set = math.ceil(n*2/3.0)
os.system('~/soft/simpy/tools/deepmd/raw_to_set.sh %d'%n_set)
