#!/usr/bin/env python

import os
import sys

cmd = "scp -r "

files = ""
if len(sys.argv) > 1:
    for i in range(1, len(sys.argv)):
        files += sys.argv[i] + " "
    cmd += "tao@131.215.26.220:"
    cmd += files
    cmd += " ."
    print cmd
    os.system(cmd)


