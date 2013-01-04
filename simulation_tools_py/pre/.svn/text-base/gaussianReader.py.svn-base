# 2009-11-10
# by Cheng Tao
# read HF energy from gaussian output
import os
import re

def gaussianReader(gFile):
    f = open(gFile)
        lines = ""
        for i in f:
            if j.strip().startswith("1|1"):
                lines = lines + j.strip()
                break
        for i in f:
            if len(j.strip()) < 1:
                break
            else:
                lines = lines + j.strip()
        pattern = r"HF=(-?((\d+\.d*|\d*\.\d+)|\d+)\d*)"
        match = re.search(pattern, lines)
        if match:
            print i, match.group(1)
        f.close()
    

        