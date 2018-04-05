"""
generate water solvated interface
"""
equil_mdp = """
title           = DFF generated gromacs input file
cpp             = /usr/bin/cpp
integrator      = md
dt              = 0.0010 ; ps !
nsteps          = 20000 ;
nstcomm         = 1
nstxout         = 2000 ; collect data every 1.0 ps
nstvout         = 0
nstfout         = 0
nstlist         = 5
cutoff-scheme   = Verlet
ns_type         = grid
rlist           = 0.5
coulombtype     = PME
rcoulomb        = 0.5
rvdw            = 0.5
DispCorr        = EnerPres
fourierspacing  = 0.12
fourier_nx      = 0
fourier_ny      = 0
fourier_nz      = 0
pme_order       = 4
ewald_rtol      = 1e-5
optimize_fft    = yes
; Berendsen temperature coupling is on in two groups
Tcoupl          = nose-hoover
tc_grps         = SOL SUB
tau_t           = 0.2 0.2
ref_t           = 300.0  0.0
; Pressure coupling is on
Pcoupl          = no ;parrinello-rahman
pcoupltype      = isotropic
tau_p           = 0.5
compressibility = 4.5e-5
ref_p           = 1.0
; Generate velocites is on at 300 K.
gen_vel         = no
gen_temp        = 300.0
gen_seed        = 94823

freezegrps          = SUB
freezedim           = Y Y Y

"""

min_mdp = """
title               =  amylose in dmso
cpp                 =  /usr/bin/cpp
;define              =  -DFLEXIBLE
constraints         =  none
integrator          =  steep
dt                  =  0.001    ; ps !
nsteps              =  1000
nstlist             =  10
ns_type             =  grid
rlist               =  0.5
coulombtype         =  Cut-off
rcoulomb            =  0.5
rvdw                =  0.5
fourierspacing      =  0.12
fourier_nx          =  0
fourier_ny          =  0
fourier_nz          =  0
pme_order           =  4
ewald_rtol          =  1e-5
optimize_fft        =  no
;
;       Energy minimizing stuff
;
emtol               =  10.0
emstep              =  0.01

freezegrps          = SUB
freezedim           = Y Y Y

"""

water = """TITLE     S  C  A  M  O  R  G
MODEL        1
ATOM      1 O    MOL     1      -0.247  -0.247  -0.247  1.00  0.00           O
ATOM      2 H    MOL     1      -0.617   0.493   0.493  1.00  0.00           H
ATOM      3 H    MOL     1       0.863  -0.247  -0.247  1.00  0.00           H
TER
ENDMDL
"""
itp = """
[ defaults ]
; nbfunc   comb-rule   gen-pairs   fudgeLJ   fudgeQQ
       1           2         yes    0.5000    0.8333
[ atomtypes ]
; name        mass      charge   ptype         sigma           epsilon
 Cu         Cu  29     63.5460       0.0000      A      0.3114       0.0209
 C          C   6      12.0107       0.0000      A      0.3431       0.4396
 H          H   1       1.0079       0.0000      A      0.2571       0.1842
 O          O   8      15.9994       0.0000      A      0.3118       0.2512
 Na         Na  11     22.9900       0.0000      A      0.2658       0.1256
 opls_116   O   8      15.9994      -0.820       A    3.16557e-01  6.50194e-01
 opls_117   H   1       1.0080       0.410       A    0.00000e+00  0.00000e+00
 wall       W   999   999.0000       0.0         A      0.20         0.10


[ bondtypes ]
;   i    j func        b0            kb

[ angletypes ]
;   i    j    k func       th0           cth

[ dihedraltypes ]
;   i    j    k    l func      phi0            cp      mult
"""

top = """
#include "uff.itp"

[ moleculetype ] 
; Name		nrexcl 
C		3

[ atoms ]
;   nr      type   resnr residue    atom    cgnr      charge        mass
     1      C       1          SUB    C      1       0.0000     12.01070

[ moleculetype ] 
; Name		nrexcl 
H		3

[ atoms ]
;   nr      type   resnr residue    atom    cgnr      charge        mass
     1      H       1          SUB    H      1       0.0000      1.00790

[ moleculetype ] 
; Name		nrexcl 
O		3

[ atoms ]
;   nr      type   resnr residue    atom    cgnr      charge        mass
     1      O       1          SUB    O      1       0.0000     15.99940

[ moleculetype ] 
; Name		nrexcl 
Cu		3

[ atoms ]
;   nr      type   resnr residue    atom    cgnr      charge        mass
     1      Cu       1          SUB    Cu      1      -0.0200     63.54600

[ moleculetype ] 
; Name		nrexcl 
Na		3

[ atoms ]
;   nr      type   resnr residue    atom    cgnr      charge        mass
     1      Na       1          SUB    Na      1       0.9600     22.99000

[ moleculetype ]
; molname   nrexcl
water     2

[ atoms ]
;   nr      type   resnr residue    atom    cgnr      charge        mass
     1  opls_116   1    SOL     O       1      -0.8476  15.9994
     2  opls_117   1    SOL     H       1       0.4238  1.00800
     3  opls_117   1    SOL     H       1       0.4238  1.00800

[ bonds ]
; i j   funct   length  force.c.
1   2   1   0.1 345000  0.1     345000
1   3   1   0.1 345000  0.1     345000


[ angles ]
; i j   k   funct   angle   force.c.
2   1   3   1   109.47  383 109.47  383

[ settles ]
; OW    funct   doh dhh
1   1   0.1 0.16330

[ exclusions ]
1   2   3
2   1   3

[ system ] 
;name 
Mol

[ molecules ]
;      Compound     #mols 
"""
import os, shutil

def add_na(infile, outfile, coords):
    f = open(infile, "r")
    lines = f.readlines()
    f.close()
    natoms = int(lines[1])
    pbc = lines[-1]
    natoms += 1
    na_line = "%5dLIG     Na%5d%8.3f%8.3f%8.3f\n"%(natoms, natoms, 
                coords[0], coords[1], coords[2])
    lines[1] = "%8d\n"%natoms
    lines[-1] = na_line
    lines.append(pbc)
    f.close()
    o = open(outfile, "w")
    for i in lines:
        o.write(i)
    o.close()

def gen_top_itp():
    f = open("POSCAR", "r")
    lines = f.readlines()
    ele = lines[5].strip().split()
    natoms = lines[6].strip().split()
    o = open("topol.top", "w")
    o.write(top)
    for i in range(len(ele)):
        o.write("%s %s\n"%(ele[i], natoms[i]))
    o.write("Na 1\n")
    o.write("water 48\n")
    o.close()
    o = open("uff.itp", "w")
    o.write(itp)
    o.close()

def gen_water_pdb():
    o = open("water.pdb", "w")
    o.write(water)
    o.close()

def gen_mdp():
    o = open("min.mdp", "w")
    o.write(min_mdp)
    o.close()
    o = open("equil.mdp", "w")
    o.write(equil_mdp)
    o.close()
    
def get_restrain():
    fixed_atoms = []
    n = 1
    f = open("POSCAR", "r")
    for i in f:
        if "Selective" in i:
            break
    for i in f:
        if "Direct" in i:
            break
    for i in f:
        if "F" in i:
            fixed_atoms.append(n)
        n += 1
    f.close()
    return fixed_atoms

os.popen("python ~/Soft/simpy/lib/e_2_contcar.py POSCAR")
os.popen("editconf -f sim.pdb -o a0.gro")
if 1:
    os.popen("editconf -f a0.gro -o a1.gro -rotate 180 0 0")
os.popen("editconf -f a1.gro -o a2.gro -box 1.0225 1.0225 2.0")
os.popen("editconf -f a2.gro -o a3.gro -translate 0 0 -0.7")
coords = [0.500, 0.500, 1.498]
add_na("a3.gro", "a4.gro", coords)
gen_water_pdb()
os.popen("genbox -cp a4.gro -ci water.pdb -nmol 48 -try 1000 -o a5.gro")
os.popen("editconf -f a5.gro -o a6.gro -box 1.02250   1.02250 4 -noc")
shutil.copy("a6.gro", "conf.gro")
fixed_atoms = get_restrain()
gen_top_itp()
gen_mdp()
os.popen("grompp -f min -maxwarn 1")
os.popen("mdrun && cp confout.gro conf.gro")
os.popen("grompp -f equil -maxwarn 1")
os.popen("mdrun && cp confout.gro conf.gro")
os.popen("editconf -f conf.gro -o a7.pdb")
shutil.copy("POSCAR", "POSCAR-vac")
os.popen("python ~/Soft/simpy/lib/e_2_pdb.py a7.pdb -c")
shutil.copy("POSCAR", "POSCAR-na-water.vasp")
cmd = "python /home/tao/Soft/simpy/tools/vasp/fix_atoms "
for i in fixed_atoms:
    cmd += "%d "%i
os.popen(cmd)
#os.rename("POSCAR", "POSCAR-tmp")


