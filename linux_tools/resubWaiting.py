#!/usr/bin/env python

import os
runId = []
runFolder = []
waitId = []
waitFolder = []

f = os.popen("job.py")
for i in f:
    if i.startswith("-------------------------------Job-Running"):
        break
for i in f:
    if i.startswith("-------------------------------Job-Waiting"):
        break
    else:
        if len(i.split()) > 1 :
            runId.append(i.split()[0])
            runFolder.append(i.strip().split()[1])
            
for i in f:
    if len(i.split()) > 1:
        waitId.append(i.split()[0])
        waitFolder.append(i.strip().split()[1])
f.close()

o = open("killRun.sh", 'w')
for i in runId:
    o.write("qdel %s\n"%i)
o.close()

o = open("resubRun.sh", 'w')
for i in runFolder:
    o.write("cd %s\n"%i)
    o.write("qsub rerun.sh\n")
o.close()


o = open("killWait.sh", 'w')
for i in waitId:
    o.write("qdel %s\n"%i)
o.close()

o = open("resubWait.sh", 'w')
for i in waitFolder:
    o.write("cd %s\n"%i)
    o.write("qsub grompp.sh\n")
o.close()
