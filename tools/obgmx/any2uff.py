import os, time, sys, shutil
from ase.io import read, write

control_em_steep = """integrator      = steep
nsteps          = 1000
nstcomm         = 1 
nstxout         = 1000
nstvout         = 0
nstfout         = 0
nstlist         = 5
ns_type         = grid
coulombtype     = cut-off ;PME
rlist           = 0.6
rcoulomb        = 0.6
rvdw            = 0.6
DispCorr        = EnerPres
fourierspacing  = 0.12
fourier_nx      = 0
fourier_ny      = 0
fourier_nz      = 0
pme_order       = 4
ewald_rtol      = 1e-5
optimize_fft    = yes
"""
control_em_cg = """integrator      = cg
nsteps          = 1000
nstcomm         = 1 
nstxout         = 1000
nstvout         = 0
nstfout         = 0
nstlist         = 5
ns_type         = grid
coulombtype     = cut-off ;PME
rlist           = 0.6
rcoulomb        = 0.6
rvdw            = 0.6
DispCorr        = EnerPres
fourierspacing  = 0.12
fourier_nx      = 0
fourier_ny      = 0
fourier_nz      = 0
pme_order       = 4
ewald_rtol      = 1e-5
optimize_fft    = yes
"""

control_nvt = """integrator      = md
dt              = 0.0010  
nsteps          = 100000
nstcomm         = 1 
nstxout         = 1000
nstvout         = 0
nstfout         = 0
nstlist         = 5
ns_type         = grid
coulombtype     = cut-off ;PME
rlist           = 0.6
rcoulomb        = 0.6
rvdw            = 0.6
DispCorr        = EnerPres
fourierspacing  = 0.12
fourier_nx      = 0
fourier_ny      = 0
fourier_nz      = 0
pme_order       = 4
ewald_rtol      = 1e-5
optimize_fft    = yes
; Berendsen temperature coupling is on in two groups
Tcoupl          = v-rescale
tc_grps         = system
tau_t           = 0.1 
ref_t           = 298.0 
; Generate velocites is on at 300 K.
gen_vel         = yes
gen_temp        = 298.0
gen_seed        = 94823
"""

def read_itp(itpfile):
    parms = {}
    term =''
    f = open(itpfile, 'r')
    for i in f:
        if i.startswith('#') or i.startswith(';'):
            pass
        elif len(i.strip()) < 1:
            pass
        elif i.startswith('['):
            term = i.split()[1]
            if not term in parms:
                parms[term] = []
        else:
            if term !='':
                clear_comments = 1
                if clear_comments:
                    if ';' in i:
                        parms[term].append(i.split(';')[0].split())
                    else:
                        parms[term].append(i.split())
    f.close()
    return parms

def write_itp(itpfile, parms):
    o = open(itpfile, 'w')
    ISOTIMEFORMAT='%Y-%m-%d %X'
    currentTime = time.strftime( ISOTIMEFORMAT, time.localtime( time.time() ) )
    #o.write("#CREATED AT%s\n\n"%currentTime)
    
    sections = ['moleculetype', 'atoms', 'bonds', 'angles', 'dihedrals']
    for i in sections:
        if i in parms.keys():
            o.write("[ %s ]\n"%i)
            for ii in parms[i]:
                line = ''
                for iii in ii:
                    line += iii + ' '
                line += '\n'
                o.write(line)
    o.close()

def get_uff(fname):
    os.system('obgmx -H -G 7 %s > uff.log'%fname)
    time.sleep(0.1)

def update_charge():
    shutil.copy('obgmx.itp', 'obgmx.itp.0')
    charges = []
    f = open('q.dat', 'r')
    for i in f:
        tokens = i.strip().split()
        if len(tokens) > 0:
            q = float(tokens[0])
            charges.append(q)
    f.close()

    itp = read_itp('obgmx.itp.0')
    if 'atoms' in itp.keys():
        for i in range(len(itp['atoms'])):
            itp['atoms'][i][6] = '%.6f'%charges[i]
    print(itp['atoms'])
    write_itp('obgmx.itp', itp)

def run_gromacs(fname):
    folder = '01-run-gromacs'
    if not os.path.exists(folder):
        os.mkdir(folder)
    shutil.copy(fname, folder)
    shutil.copy('obgmx.top', folder)
    shutil.copy('obgmx.itp', folder)
    os.chdir(folder)

    o = open('emsteep.mdp', 'w')
    o.write(control_em_steep)
    o.close()

    o = open('emcg.mdp', 'w')
    o.write(control_em_steep)
    o.close()

    o = open('equil.mdp', 'w')
    o.write(control_nvt)
    o.close()

    os.system('gmx grompp -f emsteep.mdp -c %s -p obgmx.top -maxwarn 1'%fname)
    os.system('gmx mdrun')

    os.system('gmx grompp -f emcg.mdp -c confout.gro -p obgmx.top')
    os.system('gmx mdrun')

    os.system('gmx grompp -f equil.mdp -c confout.gro -p obgmx.top')
    os.system('gmx mdrun')

    os.chdir('..')

if len(sys.argv) > 1:
    fname = sys.argv[1]

    get_uff(fname)
    if os.path.exists('q.dat'):
        update_charge()
    run_gromacs(fname)
    
