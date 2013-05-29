import os

def parsefolder(folder):
    for i in os.listdir(folder):
        fullname = os.path.join(folder, i)
        if os.path.isdir(fullname):
            parsefolder(fullname)
        elif i == "traj.trr":
            shfile = os.path.join(folder, "rdf.sh")
            o = open(shfile, 'w')
            o.write("""#!/bin/bash
#PBS -l nodes=1:ppn=1
cd %s
echo 4 3 | g_rdf -f traj.trr -s equil.tpr -n index.ndx -o sna
echo 4 5 | g_rdf -f traj.trr -s equil.tpr -n index.ndx -o sow"""%folder)
            o.close()
            os.system("qsub %s"%shfile)
if __name__ == "__main__":
    parsefolder('/home/chengtao/forpaper/forcefield/trr')
            
