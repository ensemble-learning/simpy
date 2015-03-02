coords = {
2059:[],
2245:[],
3699:[],
3997:[],
3998:[],
}

f = open("dump.lammpstrj", "r")

for i in f:
    tokens = i.strip().split()
    if len(tokens) == 8:
        id = int(tokens[0])
        x = float(tokens[2])
        y = float(tokens[3])
        z = float(tokens[4])
        if id in coords.keys():
            coords[id] = [x, y, z]
f.close()

for i in coords.keys().sort():
    print i, coords[i][0], coords[i][1], coords[i][2]
        
