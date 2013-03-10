"""Convert the "molecular" style data file to "chagre" style data file.
Note: Usually the dff generate "molecular" style. But in most of the case, 
"charge" style is enough.
"""

import shutil

def mol2charge(datafile):
    shutil.copy(datafile, datafile + ".bak")
    lines = ""
    f = open(datafile, "r")
    for i in f:
        lines.append(i)
        if i.strip().startswith("Atoms"):
            break
    for i in f:
        tokens = i.strip().split()
        if len(tokens) > 0:
            tokens[1] = ''
            line = " ".join(tokens) + "\n"
        else:
            line = i
        lines.append(line)
    f.close()
    o = open(datafile, "w")
    for i in lines:
        o.write(i)
    o.close()

mol2charge("data_lammps")
