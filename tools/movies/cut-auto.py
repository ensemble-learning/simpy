import time, os, shutil

def get_sec(time_str):
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + float(s)

def cut(filename):
    infile = filename
    outfile = infile.split(".")[0]
    cmd = "ffmpeg -i %s 2>&1 | grep Duration | cut -d ' ' -f 4 | sed s/,//"%infile
    dt = 60*5 # 5 min
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
        folder = "%s-%03d"%(outfile, i)
        if not os.path.exists(folder):
            os.mkdir(folder)
        cmd = 'avconv -i %s -ss %s -t %s  -vcodec copy -acodec copy -metadata track="1" "%s-%03d%s"'%(infile, ts0, ts1, outfile, i, file_ext)
        os.system(cmd)
        shutil.move("%s-%03d%s"%(outfile, i, file_ext), folder)
        t0 = t1

def main():
    file_ext = ".avi"
    for i in os.listdir("."):
        if i.endswith(file_ext):
            cut(i)

if __name__ == "__main__":
    main()
    
