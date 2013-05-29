import os
import os.path
dirList = os.listdir("./")

for i in dirList:
    if os.path.isdir(i):
        os.chdir(i)
        os.system("grompp_mpi -np 6 -f run.mdp -c run -p 2,3-dimethylbutane -o run")
        os.system("/opt/openmpi/bin/mpirun -np 6 mdrun_mpi -s run -c runi")
        os.chdir("../")


