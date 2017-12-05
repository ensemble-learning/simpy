#!/usr/bin/env python

import os

run_jobs = []
wait_jobs = []
total_jobs = ""

f = os.popen("showq -u")
for i in f:
    if i.strip().startswith("ACTIVE JOBS"):
        break
for i in f:
    if i.strip().startswith("WAITING JOBS"):
        break
    if i.strip().startswith("Total Jobs"):
	totoal_jobs = i
        break
    run_jobs.append(i)
for i in f:
    if i.strip().startswith("COMPLETING/ERRORED JOBS"):
        break
    if i.strip().startswith("Total Jobs"):
	totoal_jobs = i
        break
    wait_jobs.append(i)

#print run_jobs
#print wait_jobs
f.close()

if len(run_jobs) > 3:
    print "-"*20,
    print "Jobs running",
    print "-"*20
    for i in run_jobs[2:-1]:
        tokens = i.strip().split()
        jobid = tokens[0]
        print jobid, 
        f = os.popen("scontrol show job %s"%jobid)
        for j in f:
            if j.strip().startswith("RunTime"):
                tokens = j.strip().split()
                runtime = tokens[0].split("=")[1]
                break
        print runtime,
        for j in f:
            if j.strip().startswith("WorkDir"):
                tokens = j.strip().split()
                workfolder = tokens[0].split("=")[1][21:]
                workfolder = "~" + workfolder
                break
        print workfolder
        f.close()

if len(wait_jobs) > 3:
    print "-"*20,
    print "Jobs waiting",
    print "-"*20
    for i in wait_jobs[2:-1]:
        tokens = i.strip().split()
        jobid = tokens[0]
        print jobid, 
        f = os.popen("scontrol show job %s"%jobid)
        for j in f:
            if j.strip().startswith("WorkDir"):
                tokens = j.strip().split()
                workfolder = tokens[0].split("=")[1][21:]
                workfolder = "~" + workfolder
                break
        print workfolder
        f.close()

print "-"*20,
print "Jobs summary",
print "-"*20
print totoal_jobs

