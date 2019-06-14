f = open("OUTCAR", "r")

flag = 1
while(flag):
    flag = 0
    for i in f:
        flag = 1
        if i.strip().startswith("E-fermi"):
            break
    for i in f:
        if i.strip().startswith("band"):
            break
    for i in f:
        tokens = i.strip().split()
        if len(tokens) == 3:
            if float(tokens[2]) > 0:
                e_homo = float(tokens[1])
            else:
                e_lumo = float(tokens[1])
                print(e_lumo - e_homo)
                break
        else:
            break
