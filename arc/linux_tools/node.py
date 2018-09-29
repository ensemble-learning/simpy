#!/usr/bin/env python 
import os

os.system("pbsnodes > all_nodes.log")

f = open("all_nodes.log" , 'r')

cpu_remain = 1
for i in f:
    if i[:7] == "compute":
        print i[:-1]
    if i[13:-1] == "job-exclusive":
        print "no cpu"
    if i[5:9] == "jobs":
        if len(i.split(",")):
            no_re = 8 - len(i.split(","))
        else:
            no_re = 8
        print no_re
      
f.close()

os.system("rm all_nodes.log")
