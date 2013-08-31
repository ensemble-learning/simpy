""" Matripulate the reaxFF data base.
@log:
@note:
@todo:
"""
import os, shutil
import ConfigParser
import argparse

PATH_DB = "/net/hulk/home6/chengtao/soft/simpy/database/reaxFF"

class Log():
    def __init__(self,):
        self.name = ''
        self.path = ''
        self.system = ''
        self.element = []
        self.ff_now = ''

def parseConfig(log):
    Config = ConfigParser.ConfigParser()
    Config.read("config.ini")
    sections = Config.sections()
    if "Description" in sections:
        o = Config.options("Description")
        if "system" in o:
            log.system = Config.get("Description", "system") 
        if "element" in o:
            tokens = Config.get("Description", "element")
            log.element = tokens.split()
    if "Update" in sections:
        o = Config.options("Update")
        if "version" in o:
            log.ff_now = Config.get("Update", "version") 

def parseFolders():
    ffs = {}
    current = os.getcwd()
    os.chdir(PATH_DB)
    for i in os.listdir("."):
        if os.path.isdir(i):
            log = Log()
            os.chdir(i)
            folder = os.path.join(PATH_DB, i)
            log.name = i
            log.path = folder
            parseConfig(log)
            ffs[i] = log
            os.chdir("..")
    os.chdir(current)
    return ffs

def printStatus(ffs):
    print ""
    print "-"*50
    print "    Current status of the ReaxFF database"
    print "-"*50
    print ""
    print "%-8s%-40s%-40s"%("NAME", "DESCR", "ELEMENTS")
    print "%-8s%-40s%-40s"%("____", "_____", "________")
    
    for i in ffs.keys():
        print "%-8s%-40s%-40s"%(ffs[i].name, ffs[i].system, 
              " ".join(ffs[i].element))
    print ''
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-status", action="store_true", help="Print current status of data base.")
    parser.add_argument("-checkout", nargs=1, help="Check out force field parameters to current folder")
    args = parser.parse_args()
    
    ffs = parseFolders()

    if args.status:
        printStatus(ffs)
    if args.checkout:
        key = args.checkout[0]
        if key in ffs.keys():
            target = os.path.join(PATH_DB, key)
            target = os.path.join(target, ffs[key].ff_now)
            if os.path.exists("ffield"):
                print "Are you sure to cover current ffield? "
                flag = raw_input("Y/N: ")
                if flag == "Y" or flag == "y":
                    shutil.copy(target, "ffield")

if __name__ == "__main__":
    main()
