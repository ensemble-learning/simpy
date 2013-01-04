def spc_tip4p(coords):
    R1 = 0.1346
    R2 = 0.7308
    #R1 = 0.128
    #R2 = 0.744
    coords_tip4p = []
    x,y,z = [],[],[]
    
    counter = 0
    natom = 0
    for i in range(len(coords)):
        natom +=1
        if counter % 3 == 0:
            if len(x) == 0:
                pass
            else:
                xadd = R2*float(x[0]) + R1*(float(x[1])+float(x[2]))
                yadd = R2*float(y[0]) + R1*(float(y[1])+float(y[2]))
                zadd = R2*float(z[0]) + R1*(float(z[1])+float(z[2]))
                atom_add = ['ATOM', natom, 'M', 'MOL', mol, xadd, yadd, zadd, '1.00', '0.00', 'M']
                coords_tip4p.append(atom_add)
                x,y,z = [],[],[]
                natom += 1
        #ATOM      1  O   MOL     1     2.251   6.437   6.866
        x.append(coords[i].split()[5])
        y.append(coords[i].split()[6])
        z.append(coords[i].split()[7])
        coords_tip4p.append(coords[i].strip().split())
        coords_tip4p[natom-1][1] = natom
        mol = coords[i].split()[4]
        counter += 1
    
    natom += 1
    xadd = R2*float(x[0]) + R1*(float(x[1])+float(x[2]))
    yadd = R2*float(y[0]) + R1*(float(y[1])+float(y[2]))
    zadd = R2*float(z[0]) + R1*(float(z[1])+float(z[2]))
    atom_add = ['ATOM', natom, 'M', 'MOL', mol, xadd, yadd, zadd, '1.00', '0.00', 'M']
    coords_tip4p.append(atom_add)
    x,y,z = [],[],[]
    
    
    return coords_tip4p