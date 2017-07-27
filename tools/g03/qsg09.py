import argparse
import shutil

def update_gjf(ppn):
    f = open("g09.gjf", "r")
    lines = f.readlines()
    f.close()
    shutil.copy("g09.gjf", "g09-s.gjf")
    o = open("g09.gjf", "w")
    o.write(r"%nprocshared=")
    o.write("%d\n"%ppn)
    for i in lines:
        o.write(i)
    o.close()

PBS = """#!/bin/csh    
#PBS -l nodes=1:ppn=%ppn%
#PBS -l walltime=48:00:00
#PBS -l pvmem=2gb
#PBS -j oe

setenv g09root "/net/hulk/PMD/tcheng/00_soft/gaussian/09-e1"
source $g09root/g09/bsd/g09.login

setenv GAUSS_SCRDIR "/net/hulk/PMD/tcheng/tmp"
set g09 = /net/hulk/PMD/tcheng/00_soft/gaussian/09-e1/g09/g09

if ($?GAUSS_EXEDIR) then
  setenv GAUSS_EXEDIR ${GAUSS_EXEDIR}:${g09:h}
else
  setenv GAUSS_EXEDIR $g09:h
endif

if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH ${g09:h}:${LD_LIBRARY_PATH}
else
  setenv LD_LIBRARY_PATH $g09:h
endif

cd $PBS_O_WORKDIR
g09 g09.gjf
"""

ppn = 1

parser = argparse.ArgumentParser()

parser.add_argument("gjf", nargs=1, help="gjf file")
parser.add_argument("-ncpu", nargs=1, type=int, help="number of cpus")
parser.add_argument("-hours", nargs=1, type=int, help="time to run the job (h).")
args = parser.parse_args()

gjf = args.gjf[0]
if gjf != "g09.gjf":
    shutil.copy(gjf, "g09.gjf")

if args.ncpu:
    ppn = args.ncpu[0]

if ppn > 1:
    update_gjf(ppn)

pbs = PBS
pbs = pbs.replace("%ppn%", "%d"%ppn)
o = open("pbs", "w")
o.write(pbs)
o.close()

