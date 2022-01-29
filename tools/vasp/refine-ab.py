f = open('ML_AB.0', 'r')
lines = f.readlines()
f.close()

nc = 0
for n, i in enumerate(lines):
    if i.strip().startswith('Configuration num'):
        nc += 1
        lines[n] = '     Configuration num.%7d\n'%nc
        print(lines[n])
lines[4] = '%11d\n'%nc

o = open('ML_AB', 'w')
for i in lines:
    o.write(i)
o.close()


