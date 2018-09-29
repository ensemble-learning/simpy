dist = {}
f = open("molid.out", "r")

for i in f:
    tokens = i.strip().split()
    if i.strip().startswith("#"):
        pass
    else:
        if len(tokens) > 2:
            molname = tokens[2]
            if molname in dist.keys():
                dist[molname] += 1
            else:
                dist[molname] = 1
f.close()

for i in dist.keys():
    print i, dist[i]
