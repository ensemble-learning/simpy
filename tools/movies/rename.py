import os
for i in os.listdir("."):
    if i.endswith(".mp4"):
        tokens = i.strip().split()
        new_name = "%02d"%(int(tokens[0]))+"-"+"%02d"%(int(tokens[2]))
        new_name += ".mp4"
        os.rename(i, new_name)

