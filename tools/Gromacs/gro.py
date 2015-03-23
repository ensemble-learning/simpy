MASS = {"O":15.999, "H": 1.008, "Cl":35.453, "Pt":195.08}
atoms = []
f = open("conf.gro", "r")

tokens = f.readline()
note = tokens.strip()

tokens = f.readline()
natoms = int(tokens.strip())

for i in range(natoms):
    tokens = f.readline()
    n = int(tokens[:5])
    res = tokens[5:10].strip()
    ele = tokens[10:15].strip()
    nmol = int(tokens[15:20])
    x = float(tokens[20:28])
    y = float(tokens[28:36])
    z = float(tokens[36:42])
    atoms.append([n, res, ele, nmol, x, y, z])

tokens = f.readline()
pbc = tokens.strip().split()

o = open("topol.top", "w")
o.write("""; Toplogy file exported from DFF
; Include forcefield parameters
#include "ff.itp"


[ moleculetype ]
; Name                  nrexcl
SIM                     3

[ atoms ]
;   nr      type   resnr residue    atom    cgnr      charge        mass
""")

counter = 0
for i in atoms:
    o.write("%8d%8s%8d%8s%8s%8d%8.4f%12.4f\n"%
            (counter + 1, i[2], i[3], i[1], i[2], counter + 1, 0.0, MASS[i[2]]))
    counter += 1

o.write("""[ bonds ]
;   ai    aj funct              c0              c1              c2              c3

[ angles ]
;   ai    aj    ak funct                c0              c1              c2              c3

[ exclusions ]

[ system ]
; Name
SIM

[ molecules ]
;      Compound     #mols
SIM 1

""")
