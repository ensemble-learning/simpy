"""
"""
def get_species(molname):
    f = open("species.out", "r")
    o = open("%s.dat"%molname, "w")

    pos = 0
    for i in f:
        tokens = i.strip().split()
        if i.startswith("#"):
            for j in range(len(tokens)):
                if tokens[j] == molname:
                    pos = j
                    break
        else:
            nstep = int(tokens[0])
            if pos > 0:
                nmol = int(tokens[pos-1])
            else:
                nmol = 0
            pos = 0
            o.write("%d\t%d\n"%(nstep, nmol))

    f.close()

get_species("N3H")
        
