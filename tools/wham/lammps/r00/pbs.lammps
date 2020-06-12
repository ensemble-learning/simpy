#!/bin/bash

#SBATCH -A wag
#SBATCH -J jobs
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -o gromacs%j.%N.out
#SBATCH --time=2:00:00   # walltime
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks)
#SBATCH --nodes=1


cd $SLURM_SUBMIT_DIR
~/soft/lammps/lammps-3Mar20/src/lmp_mpi -in in.lmp > log
