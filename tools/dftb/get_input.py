import os, shutil
import socket
import argparse

print(socket.gethostname())

class Param():
    def __init__(self,):
        self.opt = 0
        self.md = 0
        self.excited = 0
        self.d3 = 1

# HubbardDerivs
# http://www.dftb.org/parameters/download/3ob/3ob-3-1-cc/
HD = { "C": -0.1492, "H": -0.1857, "O":-0.1575, "N":-0.1535,
        "Cu":-0.2000, "S":-0.11, "Na":-0.0454, "F":-0.1623,
     "Ca":-0.0340, "Cl": -0.0697, "Mg":-0.02, "P": -0.14,
      "K": -0.0339, }
# MaxAngularMomentum
MM = { "C": "p", "H": "s", "O":"p", "N":"p", "Cu":"d", "S":"d",
      "Na": "p", "F":"p", "Ca":"p", "Cl": "p", "Mg":"p", "P": "p", 
      "K": "p",}
# CovalentRadius 
# http://periodictable.com/Properties/A/CovalentRadius.html
CR = {"C":0.77, "H":0.37, "O": 0.73, "N":0.75, "Cu":1.38, "S":1.02,
      "Na":0.71, "F":0.71, "Ca":1.76, "Cl":0.99, "P":1.07, "K":2.03,}
# HybridPolarisations
HP = {
"C": [1.382, 1.382, 1.382, 1.064, 1.064, 1.064, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 2.50],
"H": [0.386, 0.386, 0.000, 0.000, 0.000, 0.000, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 0.80],
"O": [0.560, 0.560, 0.000, 0.000, 0.000, 0.000, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 3.15],
"N": [1.030, 1.030, 1.090, 1.090, 1.090, 1.090, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 2.82],
"Cu": [1.030, 1.030, 1.090, 1.090, 1.090, 1.090, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 2.82],
"S": [1.030, 1.030, 1.090, 1.090, 1.090, 1.090, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 2.82],
"Na": [0.386, 0.386, 0.000, 0.000, 0.000, 0.000, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 0.80],
# fake values
"K": [0.386, 0.386, 0.000, 0.000, 0.000, 0.000, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 0.80],
"F": [1.030, 1.030, 1.090, 1.090, 1.090, 1.090, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 2.82],
"P": [1.030, 1.030, 1.090, 1.090, 1.090, 1.090, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 2.82],
"Ca": [1.030, 1.030, 1.090, 1.090, 1.090, 1.090, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 2.82],
"Cl": [1.030, 1.030, 1.090, 1.090, 1.090, 1.090, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 2.82],
}

class Gen():
    def __init__(self, ):
        self.natoms = 0
        self.atoms = []
        self.cell = []

def gen_inp(gen, p):
    n = 0
    for i in os.listdir("."):
        if "dftb_in.hsd" in i:
            n += 1
    if n > 0:
        shutil.copy("dftb_in.hsd", "dftb_in.hsd.%d"%n)
            
    o = open("dftb_in.hsd", "w")
    o.write("Geometry = GenFormat {\n")
    o.write('    <<<"geo_start.gen"\n')
    o.write("}\n")
    o.write("\n")
    # For opt
    if p.opt:
        o.write("Driver = ConjugateGradient {\n")
        o.write(" MovedAtoms = 1:-1\n")
        o.write(" MaxSteps = 50000\n")
        o.write(' AppendGeometries = Yes\n')
        o.write(" MaxForceComponent = 1.0E-3\n")
        o.write("#LatticeOpt = Yes\n")
        o.write("#Pressure [Pa] = 1.0E5\n")
        o.write("}\n")
        o.write("\n")
    # For MD
    if p.md:
        o.write("Driver = VelocityVerlet {\n")
        o.write("MovedAtoms = 1:-1\n")
        o.write("Steps = 1000\n")
        o.write("TimeStep [fs] = 1.0\n")
        o.write("Thermostat = Berendsen {\n")
        o.write("Temperature [K] = 200\n")
        o.write("CouplingStrength = 2\n")
        o.write("AdaptFillingTemp = Yes\n")

        o.write("}\n")
        o.write("}\n")

    o.write("Hamiltonian = DFTB {\n")
    o.write("  SCC = Yes\n")
    o.write("  SCCTolerance = 1E-4\n")
    o.write("  ThirdOrderFull = Yes\n")
    o.write("  DampXH = Yes\n")
    o.write("  DampXHExponent = 4.00\n")
    o.write("  Filling = Fermi{\n")
    o.write("    Temperature [K] = 300\n")
    o.write("  }\n")
    o.write("  HubbardDerivs {\n")
    # HbbardDerivs parameters
    for i in gen.atoms:
        ele = i
        hd = HD[ele]
        o.write("    %s = %.4f\n"%(ele, hd))
    o.write("  }\n")
    o.write("  MaxAngularMomentum = {\n")
    # MaxAngularMomentum
    for i in gen.atoms:
        ele = i
        mm = MM[ele]
        o.write('    %s = "%s"\n'%(ele, mm))
    o.write("  }\n")
    o.write("  SlaterKosterFiles = Type2FileNames {\n")
    #o.write('    Prefix = "/net/hulk/home6/chengtao/soft/dftb/3OB-CHNOPSMgZnCaKNaFClBrI/slko/"\n')
    if socket.gethostname() == "mu01":
        o.write('    Prefix = "/opt/sourcecoude/dftb/dftb/3OB-CHNOPSMgZnCaKNaFClBrI/slko/"\n')
    else:
        o.write('    Prefix = "/home/chengtao/soft/dftb/3ob-3-1/"\n')
    o.write('    Separator = "-"\n')
    o.write('    Suffix = ".skf"\n')
    o.write("  } \n")

    if 0: # Dispersion SlaterKirkwood
        o.write("  Dispersion = SlaterKirkwood {\n")
        o.write("    PolarRadiusCharge = HybridDependentPol {\n")
        for i in gen.atoms:
            ele = i
            cr = CR[ele]
            hp = HP[ele]
            o.write("      %s = {\n"%ele)
            o.write("        CovalentRadius [Angstrom] = %.2f\n"%cr)
            o.write("        HybridPolarisations [Angstrom^3,Angstrom,] = {\n")
            for j in hp:
                o.write("%.3f "%j)
            o.write("\n")
            o.write("        }\n")
            o.write("      }\n")
        o.write("    }\n")
        o.write("  }\n")
    if p.d3: # Dispersion D3 Becke-Jonson (BJ)
        o.write("  Dispersion = DftD3 {\n")
        o.write("    Damping = BeckeJohnson {\n")
        o.write("      a1 = 0.5719\n")
        o.write("      a2 = 3.6017\n")
        o.write("    }\n")
        o.write("    s6 = 1.0\n")
        o.write("    s8 = 0.5883\n")
        o.write("  }\n")

    # Kpoints
    if gen.cell:
        o.write("  KPointsAndWeights = SupercellFolding {\n")
        o.write("   1 0 0\n")
        o.write("   0 1 0\n")
        o.write("   0 0 1\n")
        o.write("   0.0 0.0 0.0\n")
        o.write("  }\n")
    o.write("}\n")
    o.write("\n")

    if p.excited:
        o.write('ExcitedState = {\n')
        o.write('  Casida = {\n')
        o.write('    NrOfExcitations = 16\n')
        o.write('    Symmetry = singlet\n')
        o.write('  }\n')
        o.write('}\n\n')

    o.write("Analysis = {\n")
    o.write("  WriteEigenvectors = Yes\n")
    o.write("}\n")
    o.write("\n")
    o.write("Options = {\n")
    o.write("  WriteDetailedXML = Yes\n")
    o.write("}\n")
    o.write("\n")
    o.write("ParserOptions = {\n")
    o.write("  ParserVersion = 3\n")
    o.write("}\n")
    o.write("\n")
    o.close()

def read_gen(gen):
    f = open("geo_start.gen", "r")
    lines = f.readlines()
    tokens = lines[0].strip().split()
    gen.natoms = int(tokens[0])
    tokens = lines[1].strip().split()
    gen.atoms = tokens
    if len(lines) == gen.natoms + 6:
        for i in range(gen.natoms+2, gen.natoms + 6):
            tokens = lines[i].strip().split()
            gen.cell.append([float(j) for j in tokens])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("runtype", default="sp", nargs="?", help="run type")
    #parser.add_argument("-q", type=int, help="charge")
    #parser.add_argument("-nspin", type=int, help="Spin multiplicity")
    #parser.add_argument("-ncpu", type=int, help="number of cpu for DFTB calculation")
    args = parser.parse_args()

    p = Param()
    gen = Gen()
    # assign parameters
    if args.runtype == 'opt':
        p.opt=1

    if args.runtype == 'tddft':
        p.excited=1

    read_gen(gen)
    gen_inp(gen, p)
