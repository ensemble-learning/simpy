#!/usr/bin/env python2

import sys
import os
import shutil
import numpy as np

if len(sys.argv) < 3:
    print("Need three parameters!")
    print("t0 t1 nt!")
else: 
    t0 = int(sys.argv[1])
    t1 = int(sys.argv[2])
    nt = int(sys.argv[3])
    if not os.path.exists("wf"):
        os.mkdir("wf")

    data = np.linspace(t0, t1, nt)
    print t0, t1, data
    for i in data:
        fname = "POSCAR_%06d"%int(i)
        shutil.copy(fname, "wf")

