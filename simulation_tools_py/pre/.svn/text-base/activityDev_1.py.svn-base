import os
import os.path
def activityDe(kb,c):
    deltaG = (float(kb[0]) - float(kb[1]))/1000
    acc = 1/(1+c*deltaG)
    return acc
    
def parsefolder(folder):
    for i in os.listdir(folder):
        fullname = os.path.join(folder, i)
        if os.path.isdir(fullname):
            parsefolder(fullname)
        elif i.endswith(".log"):
            f = open(fullname, 'r')
            kbintergral = []
            for j in f:
                kbintergral.append(j.strip().split()[1])
            if "s30" in fullname:
                acc = activityDe(kbintergral,3.02)
            elif "s80" in fullname:
                acc = activityDe(kbintergral,6.66)
            print fullname, acc
            f.close()
            
if __name__ == "__main__":
    parsefolder('.')