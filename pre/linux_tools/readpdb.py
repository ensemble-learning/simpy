#!/usr/bin/env python
import os
import os.path

nstart = 90 
nend = 108
ninterval = 2

for a in range( nstart, nend, ninterval ):
    for b in range( nstart, nend, ninterval ):
           os.chdir("./a%03d_%03d"%(a,b))
           if os.path.getsize("md0.log") > 100000 :
               os.system("editconf_mpi -f confout.gro -o aa%03d_%03d.pdb"%(a,b))
               os.system("cp aa%03d_%03d.pdb ../pdb/"%(a,b))
           os.chdir("../")

