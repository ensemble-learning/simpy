ndx = []
f = open("index.ndx", "r")
for i in f:
    tokens = i.strip()
    ndx.append(int(tokens))
f.close()

f = open("dump.lammpstrj", "r")
o = open("new.lammpstrj", "w")

counter = 0
for i in f:
    if i.strip().startswith("ITEM"):
        o.write(i)
    else:
        tokens = i.strip().split()
        if len(tokens) == 8:
            atomid = int(tokens[0])
            if atomid in ndx:
                pass
            else:
                tokens[0] = "%d"%(counter%8487 + 1)
                line = " ".join(tokens) + "\n"
                o.write(line)
                counter += 1
        else:
            o.write(i)
o.close()
f.close()
