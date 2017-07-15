#!/usr/bin/env python

import sys
import string
import os
import argparse
import socket

# This script is useful for running Jaguar in parallel.
# With SCF calculations, efficiency considerations limits the number of
# processors you should request. The number of processors recommended is as follows:
# "Divide the number of basis functions by 100 for HF or DFT or 80 for LMP2 jobs, then
#  discard any portion of this number after the decimal place. This number is the maximum
#  number of processors advised for an efficient prallel run."
# They also recommend no more than 4 for HF or DFT, and no more than 6 for LMP2.
#
# Of course, if you're requesting several structures to be calculated, then Jaguar
# can distribute the structures between the processors more efficiently so the above
# considerations are less meaningful.

#this isn't really used...
temp_dir = "/temp1/chengtao/"

#There aren't any personalizations below this line
parser = argparse.ArgumentParser()

parser.add_argument("fname", nargs=1, help="jaguar file name")
parser.add_argument("-nnodes", nargs=1, type=int, help="number of nodes")
parser.add_argument("-ncpu", nargs=1, type=int, help="number of cpus")
parser.add_argument("-hours", nargs=1, type=int, help="time to run the job (h).")


port = 1122
nodes = 1
ppn = 1
hours = 48

args = parser.parse_args()
filename = args.fname[0]

if args.ncpu:
    ppn = args.ncpu[0]

print "We're going to use port " + str(port) + "."

processors = nodes*ppn

base = filename[:-2]
file = open(base+'run','w')

#remove the last period (from the extension)
jobname = base[:-1]
#then let the job name be the last 9 characters
#we assume that the more interesting part of the filename
#is near the end, and because qstat is likely to truncate anyway,
#we'll choose what needs to be displayed.
#jobname = "j" + jobname[-9:]

print "Using temp space: " + temp_dir

file.write("#!/bin/csh    \n")
file.write("#PBS -l nodes=%d"%nodes + ":ppn=%d"%ppn + ",walltime=%d"%hours + ":00:00,pvmem=2gb\n")
#file.write("#PBS -l nodes=node-3-4" + ":ppn=" + ppn + ",walltime=" + hours + ":00:00,pvmem=4gb\n")
file.write("#PBS -N " + jobname + " \n")
#file.write("#PBS -m ea        \n")
#file.write("#PBS -j oe         \n")
file.write("\n")

if socket.gethostname() == "hive.wag.caltech.edu":
    file.write("setenv SCHRODINGER /project/exec/schrodinger_2010\n")
else:
    file.write("setenv SCHRODINGER /exec/schrodinger_2013\n")
#file.write("setenv SCHRODINGER /project/exec/schrodinger_2012\n")
file.write("if ( -d $SCHRODINGER ) set path= ( $SCHRODINGER $path )\n")
#file.write("setenv SCHRODINGER /project/exec/schrodinger_2012\n")
#file.write("if ( -d $SCHRODINGER ) set path= ( $SCHRODINGER $path )\n")

#file.write("cd " + os.getcwd() + "\n")
file.write("cd $PBS_O_WORKDIR\n")
#file.write("cp ~/jag/methane.basis ." + "\n")
file.write("setenv JAGUAR_TEMP " + temp_dir + "\n")
file.write("setenv LM_LICENSE_FILE @10.254.1.1\n")

mpi_used = "mpich-1.2.5.10-ch_p4-gcc"
if processors == 1:
#	file.write("jaguar run -VER v75207 -WAIT " + filename + "\n")
#        file.write("setenv SCHRODINGER /project/exec/schrodinger_2010\n")

#        file.write("if ( -d $SCHRODINGER ) set path= ( $SCHRODINGER $path )\n")
        file.write("jaguar run -WAIT  " + filename + "\n") #-VER v79025 
#        file.write("jaguar run  -WAIT " + filename + "\n")
else:
	file.write("setenv SCHRODINGER_NODEFILE $PBS_NODEFILE\n")
#        file.write("setenv PATH /opt/mpich-1.2.5.10-ch_p4-gcc/bin:$PATH\n")
	file.write("$SCHRODINGER/utilities/mpich start -m $PBS_NODEFILE -p "+str(port)+"\n")
	file.write("setenv MPI_USEP4SSPORT yes\n")
	file.write("setenv MPI_P4SSPORT "+str(port)+"\n")
	file.write("cat $PBS_NODEFILE > " + base + "nodes\n")
	file.write("date\n")
        file.write("jaguar run  -WAIT -PROCS " + str(processors) + " " + filename + "\n")
#        file.write("jaguar pka -WAIT -PROCS " + str(processors) + " " + filename + "\n")

#	file.write("jaguar run -VER v75207 -WAIT -PROCS " + str(processors) + " " + filename + "\n")
	file.write("date\n")
	print "Making sure that " + mpi_used + " is selected as our MPI choice."
	os.system("switcher mpi = " + mpi_used + " --force")

#from the jaguar manual:
#the port, specified by -p, must be higher than 1023 and be a 4 digit number
#if -p is not specified, then the environment variable MPI_P4SSPORT will be used
#if MPI_P4SSPORT is not set, then the default of 1234 will be used
#how can i check which ports are in use??? can different user's choices conflict?

	if int(port) <= 1023 or int(port) > 9999:
		print "The port needs to be 1023 < port <= 9999, so " + str(port) + " is not a good choice, using something else instead."
		port = 1123

handle = os.popen("echo `hostname`")
hostname = handle.readline()
hostname = hostname.rstrip()

#email $USER the output file
#file.write("csh -c \"cat " + base + "run; cat " + base + "log;\" | ssh oscar_server '/bin/mail -s \"["+hostname+"] " + base + "log\" $USER '\n\n\n")
#file.write("date\n")
file.close()

print "Submitting " + filename + " (jobname = " + jobname + ") to " + str(processors) + " processors on " + hostname + " for " + " %d hours"%hours
command = "qsub " + base +"run"
print command
#os.system(command) 
