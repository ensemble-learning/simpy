# convert gromacs format to normal format

params = []
f = open('dihedral.itp', 'r')
for i in f:
    tokens = i.strip().split()
    if len(tokens) > 11:
        params.append(tokens)
f.close()

o = open('opls-dihedral.dat', 'w')
for i in params:
    a1, a2, a3, a4 = i[0:4]
    c0, c1, c2, c3, c4, c5 = [float(j) for j in i[5:11]]
    f4 = -0.25 * c4
    f3 = -0.5  * c3
    f2 = -c2   - c4
    f1 = -2*c1 - 1.5*c3
    o.write("%-4s%-4s%-4s%-4s  opls"%(a1, a2, a3, a4))
    o.write("%10.4f%10.4f%10.4f%10.4f\n"%(f1, f2, f3, f4))
o.close()
