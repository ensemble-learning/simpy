import sys
import numpy as np

def analysis(fname, keyword):
    data = []
    f = open(fname, "r")
    o = open("%s.xvg"%keyword, "w")
    for i in f:
        if keyword in i:
            tokens = i.strip().split("=")
            for j in range(len(tokens)):
                if keyword in tokens[j]:
                    a = tokens[j+1].strip().split()[0]
                    o.write("%s\n"%a)
                    data.append(float(a))
    o.close()
    f.close()
    data = np.array(data)
    start = int(len(data) * 0.1)
    print np.average(data[start:])
    e
if __name__ == "__main__":
    
    if len(sys.argv) > 2:
        fname = sys.argv[1]
        keyword = sys.argv[2]
        analysis(fname, keyword)
    else:
        print "readLog.py logfile keywords"
    
    
    