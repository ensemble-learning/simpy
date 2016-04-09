import os, shutil
# HubbardDerivs
HD = { "C": -0.1492, "H": -0.1857, "O":-0.1575, "N":-0.1535,
        "Cu":-0.2000, "S":-0.11, "Na":-0.0454, "F":-0.1623,
    }
# MaxAngularMomentum
MM = { "C": "p", "H": "s", "O":"p", "N":"p", "Cu":"d", "S":"d",
      "Na": "p", "F":"p",}
# CovalentRadius 
# http://periodictable.com/Properties/A/CovalentRadius.html
CR = {"C":0.77, "H":0.37, "O": 0.73, "N":0.75, "Cu":1.38, "S":1.02,
      "Na":0.71, "F":0.71,}
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
"F": [1.030, 1.030, 1.090, 1.090, 1.090, 1.090, 3.8, 3.8, 3.8, 3.8, 3.8, 3.8, 2.82],
}

class Gen():
    def __init__(self, ):
        self.natoms = 0
        self.atoms = []
        self.cell = []

def gen_inp_opt(gen):
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
    o.write("Driver = ConjugateGradient {\n")
    o.write(" MovedAtoms = 1:-1\n")
    o.write(" MaxSteps = 50000\n")
    o.write(' AppendGeometries = Yes\n')
    o.write(" MaxForceComponent = 1.0E-3\n")
    o.write("}\n")
    o.write("\n")
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
    o.write('    Prefix = "/net/hulk/home6/chengtao/soft/dftb/3OB-CHNOPSMgZnCaKNaFClBrI/slko/"\n')
    o.write('    Separator = "-"\n')
    o.write('    Suffix = ".skf"\n')
    o.write("  } \n")
    o.write("  Dispersion = SlaterKirkwood {\n")
    o.write("    PolarRadiusCharge = HybridDependentPol {\n")
    # Dispersion
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
    o.write("Options = {\n")
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

def main():
    gen = Gen()
    # assign parameters
    read_gen(gen)
    gen_inp_opt(gen)

if __name__ == "__main__":
    main()
