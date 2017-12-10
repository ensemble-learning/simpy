#!/usr/bin/env python

import os
import sys
import shutil

f = open("jobinfo", "r")
lines = f.readlines()
job_path = lines[-1].strip()
os.system("cp %s . -r"%job_path)
print job_path


