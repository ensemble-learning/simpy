def outIndex(title, sequence):
    o = open("index.ndx", 'a')
    o.write("[ %s ]"%title)
    counter = 0
    for i in sequence:
        if counter%15 == 0:
            o.write("\n")
            o.write("%-5d"%i)
        else:
            o.write("%-5d"%i)
        counter += 1
    o.write("\n")

def gOrder():
    for n in range(1,15):
        t = "c%02d"%n
        s = range(5118 + n, 5894 + n, 17)
    return t, s 

def gAngle():
    ts = [] 
    ss = []
    for n in range(10):
        ts.append("angle%02d"%n)
        ss.append([])
        for i in range(5119 + n , 5895 + n , 17):
            ss[n].append(i)
            ss[n].append(i+1)
            ss[n].append(i+2)
            ss[n].append(i+3)
    return ts, ss
    
if __name__ == "__main__":
    ts, ss = gAngle()
    for i in range(len(ts)):
        outIndex(ts[i], ss[i])