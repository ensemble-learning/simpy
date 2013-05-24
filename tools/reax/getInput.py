""" This code is used to generate reaxFF lammps_input file
according to lammps.data (also generated from simpy
"""

from template import INP


#C H O N
# @note: This is only specified for current ffield
FF = {"C": 1, "H": 2, "O": 3, "N": 4, "Ca":4, "Al":6}
MASS = {12.011:"C", 14.007: "N", 15.994:"O", 1.0079:"H", 40.078:"Ca",\
        26.982:"Al"}

def parseData(fname="lammps.data"):
    """parse the data file to get the mass infor
    """
    m = []
    f = open(fname, "r")
    for i in f:
        if i.strip().startswith("Masses"):
            break
    for i in f:
        if len(i.strip()) == 0:
            pass
        elif i.strip().startswith("Atoms"):
            break
        else:
            tokens = i.strip().split()
            amass = float(tokens[1])
            m.append(amass)
    f.close()
    return m

m = parseData()
ty = []
for i in m:
    at = MASS[i]
    ff = FF[at]
    ty.append("%d"%ff)

lines = INP.replace("%ffield_atoms%", " ".join(ty))
o = open("lammps_input", "w")
o.write(lines)
o.close()
