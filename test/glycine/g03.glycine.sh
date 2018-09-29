#!/bin/bash

#PBS -q fast
#PBS -l nodes=1:ppn=4

echo PBS work @ $PBS_NODEFILE
JOB_HOME=/home/chengtao/qm/reax2/chno/glycine
JOB_DIR=/state/partition1/jobs/${USER}
mkdir ${JOB_DIR}
JOB_DIR=${JOB_DIR}/GAU
mkdir ${JOB_DIR}
JOB_DIR=${JOB_DIR}/35811
while [ -d $JOB_DIR ];do
    JOB_DIR=${JOB_DIR}X
done
mkdir ${JOB_DIR}
export GAUSS_EXEDIR=/state/partition1/apps/g03sun:/state/partition1/apps/g03sun/bsd
export GAUSS_ARCHDIR=${JOB_DIR}
export PATH=$PATH:/state/partition1/apps/g03sun
export LD_LIBRARY_PATH=/state/partition1/apps/g03sun
export GAUSS_SCRDIR=${JOB_DIR}
export G03BASIS=/state/partition1/apps/g03sun/basis

cd ${JOB_DIR}
cp $JOB_HOME/freq.gjf .
g03 freq.gjf
rm -f *.rwf
cp $JOB_DIR/*.* $JOB_HOME
cd $JOB_HOME
rm -rf $JOB_DIR
