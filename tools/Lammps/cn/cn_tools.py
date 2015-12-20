from operator import itemgetter

def get_coords():
    """
    Read ndx file.
    """
    ndx = []
    f = open("ndx.dat", "r")
    for i in f:
        tokens = i.strip().split()
        if len(tokens) > 0:
            for j in range(len(tokens)):
                if tokens[j] != "&":
                    ndx.append(int(tokens[j]))
    
    f = open("lammps.data", "r")
    for i in f:
        if i.strip().startswith("Atom"):
            break
    
    selected = []
    for i in f:
        tokens = i.strip().split()
        if len(tokens) == 6:
            id = int(tokens[0])
            x = float(tokens[3])
            y = float(tokens[4])
            z = float(tokens[5])
            if id in ndx:
                selected.append([id,x,y,z])
    
    sorted(selected, key=itemgetter(1))
    o = open("coords.dat", "w")
    o.write("#%7s%10s%10s%10s\n"%("id", "x", "y", "z"))
    for (i,j,k,l) in selected:
        o.write("%8d"%i)
        o.write("%10.4f%10.4f%10.4f\n"%(j, k, l))
    o.close()

def main():
    get_coords()

main()
