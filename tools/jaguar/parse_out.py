"""
Parse the Jaguar out file
"""
import sys

class OneStep():
    def __init__(self,):
        self.coords = []
        self.energy = 0.0

def parse_out(fname):
    f = open(fname, "r")
    for i in f:
        if i.strip().startswith("Input geometry"):
            break

    counter = 0
    frame = OneStep()
    for i in f:
        if len(i.strip()) == 0:
            break
        else:
            if counter > 1:
                frame.coords.append(i)
        counter += 1
    
    frames = []
    flag = 1
    while(flag):
        flag = 0
        for i in f:
            flag = 1
            tokens = i.strip().split()
            if i.strip().startswith("geometry optimization step"):
                frames.append(frame)
                frame = OneStep()
            elif i.strip().startswith("energy:"):
                frame.energy = float(tokens[1])
            elif i.strip().startswith("new geometry:"):
                break
        counter = 0
        for i in f:
            if len(i.strip()) == 0:
                break
            else:
                if counter > 1:
                    frame.coords.append(i)
            counter += 1
    frames.append(frame)
    return frames

def toXYZ(frames):
    counter = 0
    for i in frames:
        o = open("config_%03d.xyz"%counter, "w")
        o.write("%d\n"%len(i.coords))
        o.write("%.6f\n"%i.energy)
        for j in i.coords:
            o.write(j)
        o.close()
        counter += 1

def usage():
    sys.stdout.write("python parse_out.py outfile\n")

def main():
    if len(sys.argv) < 2:
        usage()
    else:
        fname = sys.argv[1]
        frames = parse_out(fname)
        toXYZ(frames)

if __name__ == "__main__":
    main()
