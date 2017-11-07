import time
import os

def get_sec(time_str):
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + float(s)

infile = raw_input("input file name: ")
outfile = infile.split(".")[0]
cmd = "ffmpeg -i %s 2>&1 | grep Duration | cut -d ' ' -f 4 | sed s/,//"%infile
dt = raw_input("what is the time step (in min)?: ")
dt = int(dt)*60
token = ''
for i in os.popen(cmd):
    token = i.strip()
t_total = get_sec(token)
nc = int(t_total/dt) + 1
print nc

t0 = 0
time.strftime('%H:%M:%S', time.gmtime(0))
for i in range(nc):
    t1 = t0 + dt
    ts0 = time.strftime('%H:%M:%S', time.gmtime(t0))
    ts1 = time.strftime('%H:%M:%S', time.gmtime(dt))
    cmd = 'avconv -i %s -ss %s -t %s  -vcodec copy -acodec copy -metadata track="1" "%s-%03d.mp4"'%(infile, ts0, ts1, outfile, i)
    os.system(cmd)
    t0 = t1
