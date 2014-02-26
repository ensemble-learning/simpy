#!/usr/bin/env python

"""dump the pathes
"""
import os
import shutil
from subprocess import call

def usage():
    print """mdump: 
    usage: mdump 
    """

def main():
    lines = []
    os.chdir("/net/hulk/home6/chengtao/")
    f = open(".bookmarks", "r")
    for i in f:
        tokens = i.strip().split()
        line = i
        if len(tokens) == 3:
            line = 'alias\t%s\t"cd %s"\n'%(tokens[0], tokens[2])
        lines.append(line)
    f.close()
    shutil.copy(".bookmarks", ".bookmarks_bak")
    o = open(".bookmarks", "w")
    for i in lines:
        o.write(i)
    o.close()
   
if __name__ == "__main__":
    main()
