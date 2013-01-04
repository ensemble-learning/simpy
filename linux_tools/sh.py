#!/usr/bin/env python

import os
import sys

def usage():
    print "dffjob nodes[n] msdfile[name]"

if len(sys.argv) > 2:
    nodes = int(sys.argv[1])
    msdfile = sys.argv[2]
    
    o = open("grompp.sh", 'w')

    o.write("""#!/bin/bash
#PBS -l nodes=1:ppn=%d
cd $PBS_O_WORKDIR
"""%nodes)

    o.write("grompp_mpi -np %d -f GromacsMinimize.mdp -c %s -p %s -o em\n"%(nodes, msdfile, msdfile))
    o.write("/opt/openmpi/bin/mpirun -np %d mdrun_mpi -s em -c em\n"%nodes)

    o.write("grompp_mpi -np %d -f GromacsLiquidBulkEq.mdp -c em -p %s -o equil\n"%(nodes, msdfile))
    o.write("/opt/openmpi/bin/mpirun -np %d mdrun_mpi -s equil -c equil\n"%nodes)

    o.write("grompp_mpi -np %d -f GromacsLiquidBulkRun.mdp -c equil -p %s -o run\n"%(nodes, msdfile))
    o.write("/opt/openmpi/bin/mpirun -np %d mdrun_mpi -s run -c run\n"%nodes)

    o.close()

    os.system("qsub grompp.sh")
else:
    usage()

