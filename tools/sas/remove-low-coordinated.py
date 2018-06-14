f = open("cn.dat", "r")
lines = f.readlines()
f.close()

f = open("out.xyz", "r")
coords = f.readlines()
f.close()

new_coords = []

for i in range(len(lines)):
    nc = int(lines[i])
    if nc > 5:
        new_coords.append(coords[2+i])

o = open("out.new.xyz", "w")
o.write("%d\n\n"%len(new_coords))
for i in new_coords:
    o.write(i)
o.close()

