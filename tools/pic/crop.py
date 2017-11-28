import os, sys

lx = 360
ly = 300
x0 = 600
y0 = 150

for i in sys.argv[1:]:
    fname = i
    cmd = "convert %s -crop %dx%d+%d+%d ../%s"%(fname, lx, ly, x0, y0, fname)
    os.system(cmd)

