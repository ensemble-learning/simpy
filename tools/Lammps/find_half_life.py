nhalf = 15

f = open("N3H.dat", "r")

for i in f:
    tokens = i.strip().split()
    if len(tokens) == 2:
        step = int(tokens[0])
        nmol = int(tokens[1])
        if nmol <= nhalf:
            print step
            break
    
f.close()
