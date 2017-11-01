import os
for i in os.listdir("."):
    if i.endswith(".mp4"):
        cmd = "ffmpeg -i %s 2>&1 | grep Duration | cut -d ' ' -f 4 | sed s/,//"%i
        os.system(cmd)
