#!/usr/bin/env python

import os
import os.path
import random
import time

nodes_permit = [
"compute-0-0",
"compute-0-1",
"compute-0-10",
"compute-0-2",
"compute-0-3",
"compute-0-4",
"compute-0-6",
"compute-0-7",
"compute-0-8",
"compute-0-9",
"compute-1-0",
"compute-1-1",
"compute-1-11",
"compute-1-12",
"compute-1-14",
"compute-1-16",
"compute-1-17",
"compute-1-18",
"compute-1-2",
"compute-1-3",
"compute-1-5",
"compute-1-6",
"compute-1-7",
"compute-1-8",
"compute-1-9",
]

JOB_COUNTER = 0
def nodes_free( nodes = 8 ):
# This function is used to find available nodes in permitted nodes list

    nodes_list =[]

    nodes_infor = os.popen("freenodes", 'r')
    for i in nodes_infor:
        a = i.strip().split(":")
        if int(a[-1]) >= nodes and a[0] in nodes_permit:
            nodes_list.append(a[0])

    if len(nodes_list) > 0:
        return nodes_list
    else:
        #temp = random.randint(0,(len(nodes_permit)-1))
        return nodes_permit


def sub_all(folder, freenodes_list, assign):
    #This code is used to sub all the sh file under the folder
    global JOB_COUNTER
    listall = os.listdir(folder)
    for i in listall:
        current_dir = os.path.join(folder, i)
        if os.path.isdir(current_dir):
            sub_all(current_dir, freenodes_list, assign)
        else:
            if i[-3:] == ".sh":
                os.chdir(folder)
                if assign == "yes":
                    nodesNo = freenodes_list[ JOB_COUNTER % len(freenodes_list) ]
                    os.system("qsub -l nodes=%s:ppn=%s %s"%(nodesNo, "8", i))
                    print "sub %s %s at %s"%(folder, i, nodesNo)
                    time.sleep(0.5)
                    JOB_COUNTER += 1
                elif assign == "no":
                    os.system("qsub -l nodes=1:ppn=%s %s"%("8", i))
                    time.sleep(0.8)

if __name__ == "__main__":
    root_path = os.getcwd()
    #freenodes_list = nodes_free()
    freenodes_list = ['compute-1-6', 'compute-1-1', 'compute-1-8', 'compute-1-12', 'compute-1-14', 'compute-1-5', 'compute-0-8', 'compute-0-10', 'compute-1-18' ]
    sub_all(root_path, freenodes_list, "no")
