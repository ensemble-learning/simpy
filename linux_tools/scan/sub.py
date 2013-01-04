import os

doc_dir = os.getcwd()

n = 11
step_n = 2
split_n = n/step_n
first_n = 95

iter_n = 2
core_n = 8

for i in range( split_n ):
    o = open("sub%d.sh"%i ,'w')
    o.write("#!/bin/bash\n#PBS -l nodes=1:ppn=%d\n"%core_n)
    for a in range(first_n , first_n + step_n):
        print a
        for b in range(95,106):
               o.write(doc_dir+"/a%03d_%03d\n"%(a,b))
               for j in range(iter_n):
                   o.write("grompp_mpi -np %d -f equil.mdp -c confout.gro -p ethanol -o equil\n"%core_n)
                   o.write("mpirun -np %d mdrun_mpi -s equil\n"%core_n)
    first_n = first_n + step_n 
o.close()

