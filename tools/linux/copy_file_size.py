"""
Copy file less than 10 MB
"""
import sys
import os
import shutil

t0 = "/home/tao/dropbox/"
for root, dirs, files in os.walk(t0):
    basename = os.path.basename(root)
    newroot = "." + root[len(t0):]
    if not os.path.exists(newroot):
        os.system("mkdir -p %s"%newroot)
    for file in files:
	f0 = os.path.join(root, file)
	f1 = os.path.join(newroot, file)
	f0_size = os.path.getsize(f0)
        f0_size = f0_size/1024.0/1024.0
        if f0_size < 10:
            try:
	        shutil.copy2(f0, f1)
            except IOError:
                pass

