"""
Calculation the free energy corrections using MOPAC output.
"""
import sys

CAL2KCAL = 1000.0
J2KCAL = 4184.0

f = open(sys.argv[1], "r")

Ezpe = 0.0
T = 298
RT = 8.314*298/J2KCAL

for i in f:
    if i.strip().startswith("ZERO POINT ENERGY"):
        tokens = i.strip().split()
        Ezpe = float(tokens[3])
    if i.strip().startswith("CALCULATED THERMODYNAMIC PROPERTIES"):
        break

for i in f:
    if "VIB" in i:
        tokens = i.strip().split()
        E0_v = float(tokens[3])
        S_v = float(tokens[5])
        break

for i in f:
    if i.strip().startswith("TOT"):
        tokens = i.strip().split()
        E0 = float(tokens[2])
        S = float(tokens[4])
        break

H = E0/CAL2KCAL + Ezpe + RT
G = H - T*S/CAL2KCAL
H_v = E0_v/CAL2KCAL + Ezpe + RT
G_v = H_v - T*S_v/CAL2KCAL 

print "E_zpe: zero point energy (kcal/mol) = %.2f"%Ezpe
print "H: Enthalpy (kcal/mol) = %.2f"%H
print "G:Free Energy (kcal/mol) = %.2f"%G
print "G_v:Free Energy with only vib (kcal/mol) = %.2f"%G_v


