#!/usr/bin/env python
#This code is writen to calculate the average center of mass of TATB crysl

import sys
import os

def generateXVG():
    for i in range(1,129):
        os.system("echo %d | g_traj_mpi -f traj.trr -s equil.tpr -n index -com -ox %d"%(i,i))

def readcoord(a, file_name):
    f = open( file_name, 'r')
    counter = 0
    for i in f:
        if counter > 20:
            a.append(float(i.split()[1]))
            a.append(float(i.split()[2]))
            a.append(float(i.split()[3]))
        counter += 1
def averageall(a, b):
    counter = 0
    x = 0
    y = 0
    z = 0
    for i in a:
       if counter%3 == 0:
           x += i
       if counter%3 == 1:
           y += i
       if counter%3 == 2:
           z += i
       counter += 1
    x = x / len(a) * 3
    y = y / len(a) * 3
    z = z / len(a) * 3
    b.append(x)
    b.append(y)
    b.append(z)

def outfile(b):
    counter = 0
    outfile = open("averageCom.log" ,'w')
    for i in b:
        if counter%3 == 0:
            outfile.write('\n')
        outfile.write("%15.5f"%i)
        counter += 1

    
a = []
b = []
generateXVG()
for i in range(1,129):
    a = []
    readcoord(a,"%d.xvg"%i)
    averageall(a, b )
outfile(b)
