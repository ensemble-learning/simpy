import os
import shutil

def a1000Kmdp(folder, nst):
    f = open('run.mdp', 'r')
    fullname = os.path.join(folder, 'a1000.mdp')
    o = open(fullname, 'w')
    for i in f:
        if i.strip().startswith('nsteps'):
            steps = 1000 + nst*1000
            i = "nsteps          = %d ;\n"%steps
        elif i.strip().startswith('ref_t'):
            i = "ref_t           =  1200.0\n"
        else:
            pass
        o.write(i)
    o.close()
    f.close()
    
def shFile(folder):
    fullname = os.path.join(folder, 'grompp.sh')
    o = open(fullname, 'w')
    o.write("""#!/bin/bash
cd $PBS_O_WORKDIR
/share/apps/gromacs333n/bin/grompp_mpi -f a1000.mdp -c run -p *.top -o sample
/share/apps/gromacs333n/bin/mdrun_mpi -s sample -c sample
/share/apps/gromacs333n/bin/grompp_mpi -f run.mdp -c sample -p *.top -o equil
/share/apps/gromacs333n/bin/mdrun_mpi -s equil -c equil
/share/apps/gromacs333n/bin/grompp_mpi -f run.mdp -c equil -p *.top -o run
/share/apps/gromacs333n/bin/mdrun_mpi -s run -c run
echo Potential | g_energy > potential.log
rm *.edr *.tpr *.log \#*
""")

def copyFile(folder, nst, model):
    a1000Kmdp(folder, nst)
    shFile(folder)
    shutil.copy("run.gro", folder)    
    shutil.copy("run.mdp", folder)    
    shutil.copy(model + ".top", folder)    
    shutil.copy(model + ".itp", folder)    

def parseFolder(folder):
    for i in os.listdir(folder):
        if i.endswith(".top"):
            model = i.split('.')[0]
    for i in range(10000):
        dirname = 'md_%05d'%i
        if os.path.isdir(dirname):
            pass
        else:
            os.mkdir(dirname)
        fullname = os.path.join(folder, dirname)
        copyFile(fullname, i, model)
            
if __name__ == "__main__":
    parseFolder('./')
