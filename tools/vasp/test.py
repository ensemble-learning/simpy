POT = {"N":"N", "O":"O", "H":"H", "C":"C", "Li":"Li_sv", "S":"S", "Ti":"Ti", "P":"P",
              "Ca":"Ca_pv", "Al":"Al", "Cu":"Cu", "Na":"Na_pv", "Cl":"Cl", "Cr":"Cr", "Ga":"Ga",
               "Br":"Br", "D": "H", "Si": "Si", "Ni": "Ni", "Pt":"Pt_pv", "Co":"Co", "Cr":"Cr",
              "I":"I", "K":"K_pv", "F":"F", "W":"W", "Au":"Au", "Cs":"Cs_sv", "Mg":"Mg_pv",
              "Ag":"Ag", "Se":"Se", "B":"B", "He":"He", "Ar":"Ar", "Xe":"Xe", "Kr":"Kr",
              "Mo": "Mo", "Fe": "Fe", "As": "As", "Ge": "Ge", "Sc": "Sc", "Zr": "Zr",
              "Y": "Y_sv", "Zn": "Zn", "Cd": "Cd", "V": "V", "Mn": "Mn", "Co": "Co",
              "Nb": "Nb_sv", "Tc": "Tc", "Ru": "Ru", "Rh": "Rh", "Pd": "Pd", "Ne": "Ne",
              "Bi": "Bi", "Ba": "Ba_sv", "Hf": "Hf", "In":"In", "Ir":"Ir", "Lu":"Lu", "Os":"Os",
              "Pb":"Pb", "Te":"Te", "Re":"Re", "Sb":"Sb", "Sn":"Sn", "Sr":"Sr_sv", "Ta": "Ta", "Te":"Te",
              "Tl":"Tl", "Po":"Po", "Be":"Be_sv", "Ac": "Ac"}
keys = [i for i in POT.keys()]
keys.sort()
for i in keys:
    print(i)
