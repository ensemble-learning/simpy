#!/usr/bin/env python

f = open("molfra.out", 'r')

siCompound = {}
storage = -1

for i in f:
    if len(i.split()) == 5 and i.split()[0].isdigit():
        counter = int(i.split()[0])
        if 'Si' in i.split()[3]:
            if counter == storage:
                pass
            else:
                storage = counter
                siCompound[storage] = []
            siCompound[storage].append(i.split()[-1])

f.close()

keys = siCompound.keys()
keys.sort()

o = open('siTims.csv', 'w')
for i in keys:
    o.write(str(i))
    for j in siCompound[i]:
        o.write(',' + j)
    o.write('\n')
o.close()
