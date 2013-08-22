"""constants
"""

# Atom number - Atom name map
ELEMENT = {1:"H", 6:"C", 7:"N",  8:"O", 12: "Mg", 13: "Al", 14:"Si", 20:"Ca", 22:"Ti"}
# Atom name - Atom mass map
ELEMENT2MASS = {"H": 1.0079, "O": 15.994, "N": 14.007, "Al": 26.982, "AL": 26.982, "Ca": 40.078, "CA": 40.078, 
                "Ti":47.867, "TI":47.867, "Si": 28.086, "SI": 28.086, "Mg": 24.305, "Mg": 24.305, "C":12.011}
# Atom name - Atom mass map
MASS2ELMENT = {1: "H",12: "C", 14:"N", 16: "O", 24: "Mg",  27: "Al", 28:"Si", 40:"Ca", 32: "S", 48: "Ti"}

ATOMIC_MASS_CONSTANT = 1.6605389e-27 #kg
KG2G = 1000
