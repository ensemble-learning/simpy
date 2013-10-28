"""constants
"""

# Atom number - Atom name map
ELEMENT = {1:"H", 3:"Li", 6:"C", 7:"N",  8:"O", 9:"F", 12: "Mg", 13: "Al", 14:"Si", 
            15:"P", 16:"S", 17:"Cl", 20:"Ca", 22:"Ti", 32:"Ge"}
# Atom name - Atom mass map
ELEMENT2MASS = {"H": 1.0079, "O": 15.994, "N": 14.007, "Li":6.941, "Al": 26.982, "AL": 26.982, "Ca": 40.078, "CA": 40.078, 
                "Ti":47.867, "TI":47.867, "Si": 28.086, "SI": 28.086, "Mg": 24.305, "Mg": 24.305, "C":12.011,
                "Cl":35.453, "P":30.974, "S":32.065, "Ge":72.64, "GE":72.64, "B": 10.811, "ZR":91.224, "V": 50.942, "Mo":95.94,
                "MO": 95.94, "Te": 127.6, "TE": 127.6, "Nb":92.906, "NB":92.906}
# Atom name - Atom mass map
MASS2ELMENT = {1: "H",12: "C", 14:"N", 16: "O", 24: "Mg",  27: "Al", 28:"Si", 32: "S", 40:"Ca", 48: "Ti"}

ATOMIC_MASS_CONSTANT = 1.6605389e-27 #kg
KG2G = 1000
