#!/usr/bin/env python
import os

for i in range(84, 100,4):
    for j in range(84, 100, 4):
        os.chdir("a%03d_%03d"%(i,j))
        os.system("echo %d %d >> ../cryinfor.log"%(i,j))
        os.system("cry.py >> ../cryinfor.log")
        os.chdir("../")
