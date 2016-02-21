#!/usr/bin/env python

"""
Generate pbs script for JOBS calculation on HPC.
exJOBS ncpu

@log:    1/18/2015: argparse
@log:    2/10/2016: slurm
"""

import sys
import os
import shutil
import argparse

# System dependent information
N_CPU_PER_NODE = 16
PBS_QUEUE = "workq"
EXE_LOC = "/net/hulk/home6/chengtao/soft/vasp.5.3/vasp.5.3/"
RUN_EXE = "vasp"
OUTPUT_FOLDER = 
PRE_LOAD = """
# preload lib
source /share/apps/intel/composer_xe_2013_sp1.3.174/composer_xe_2013_sp1.3.174/bin/compilervars.csh intel64
setenv OPAL_PREFIX /project/source/openmpi/build-1.8.3-ion
setenv PATH ${OPAL_PREFIX}/bin:${PATH}
setenv LD_LIBRARY_PATH ${OPAL_PREFIX}/lib:$LD_LIBRARY_PATH
"""

class Vaspjob():
    def __init__(self,):
        self.nnode = 0
        self.ncpu = 0
        self.vasploc = ""
        self.vaspexe = ""
        self.queue = ""

INP = """#PBS -V
#PBS -j oe
#PBS -N JOBS
#PBS -q %pbsqueue%
#PBS -l nodes=1:ppn=%ncpu%,walltime=24:00:00
#PBS -S /bin/tcsh

set MPI=mpirun
set EXELOC=%EXELOC%
set exe=$EXELOC/%RUNEXE%

echo "The JOBS executable is: " $exe

set DATADIR=$PBS_O_WORKDIR
set NPROCS=`wc -l $PBS_NODEFILE | cut -f1 -d " "`
set JOBID=`echo $PBS_JOBID | awk -F'.' '{print $1}'`
set OUTPUTDIR=/temp1/JOBS/$JOBID

echo "The JOBS output will temporarily be written to :" $OUTPUTDIR
echo "JOBS Job: " $JOBID "is starting with " $NPROCS "cores on:"
cat $PBS_NODEFILE
echo

# Create Directories
#
mkdir -p $OUTPUTDIR

#
# export job info
#
echo $PBS_JOBID >  $DATADIR/jobinfo
cat $PBS_NODEFILE >>  $DATADIR/jobinfo
echo $OUTPUTDIR >>  $DATADIR/jobinfo

#
# Copy executable and any required input files to your work directory
#
cp $DATADIR/* $OUTPUTDIR
cd $OUTPUTDIR

#
# Run JOBS
#
date
$MPI -np $NPROCS -machinefile $PBS_NODEFILE $exe > LOG
echo
date
echo "JOBS has finished."
echo "The JOBS output is being copied back to: " $PBS_O_WORKDIR
cd ..
tar -cvzf $JOBID.tar.gz $JOBID
cp -r $JOBID.tar.gz $PBS_O_WORKDIR
rm $JOBID.tar.gz

exit 0
"""

def update_npar(runjob):
    ncpu = runjob.ncpu

    if ncpu < 512 and ncpu > 16 :
        npar = ncpu/16
    elif ncpu <= 16:
        npar = 2
    else:
        npar = 16

    lines = []
    f = open("INCAR", "r")
    for i in f:
        line = i
        if "NPAR" in i:
            line = " NPAR = %d\n"%npar
        lines.append(line)
    f.close()
    shutil.copy("INCAR", "INCAR.bak")
    o = open("INCAR", "w")
    for i in lines:
        o.write(i)
    o.close()

def write_input(runjob):

    o = open("runJOBS.in", "w")
    lines = INP
    lines = lines.replace("%pbsqueue%", "%s"%runjob.queue)
    lines = lines.replace("%nnode%", "%d"%runjob.nnode)
    lines = lines.replace("%ncpu%", "%d"%runjob.ncpu)
    lines = lines.replace("%EXELOC%", "%s"%runjob.vasploc)
    lines = lines.replace("%EXE%", "%s"%runjob.vaspexe)
    o.write(lines)
    o.close()
    os.system("chmod +x runJOBS.in")

def main(args):
    runjob = Vaspjob()

    runjob.ncpu = N_CPU_PER_NODE
    runjob.queue = PBS_QUEUE

    if args.ncpu:
        runjob.ncpu = int(args.ncpu[0])
    runjob.nnode = runjob.ncpu/N_CPU_PER_NODE

    exeloc = EXE_LOC
    runexe = "vasp"

    if args.runexe:
        n = args.runexe[0]
        if n == 1:
            runexe = "vasp.full"
        elif n == 2:
            runexe = "vasp.gamma"
    sys.stdout.write("Using version of JOBS: %s"%exeloc)
    sys.stdout.write("%s\n"%runexe)
    sys.stdout.flush()

    runjob.exeloc = exeloc
    runjob.runexe = runexe
        
    write_input(runjob)
    # calculate NPAR for parallel calculation
    #update_npar(runjob)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-ncpu", nargs=1, type=int, help="number of cpus")
    parser.add_argument("-runexe", nargs=1, type=int, help="vasp version: 1: v53 full, 2: v53 gamma only")
    args = parser.parse_args()
    main(args)
