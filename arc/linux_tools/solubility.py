import os
def runGromacs():
    for i in os.listdir("./"):
        if i.startswith("a") and i[-3:] == "gro":
            os.system("grompp_mpi -f equil.mdp -c %s -p *.top -o single"%i)
            os.system("mdrun_mpi -s single -g md_%04d"%int(i[1:5]))
            os.system("rm \#*")

def genGro():
    f1 = open("coord.xvg", 'r')
    line_no = 0
    for i in f1:
        if len(i.split()) == 1507:
            a = i.split()
            f2 = open("template.gro", "r")
            o = open("a%04d.gro"%line_no, 'w')
            counter = 0
            for j in f2:
                if len(j.split()) == 9:
                    line = j[:20] + "%8.3f%8.3f%8.3f\n"%(float(a[counter + 1]), float(a[counter + 2]), float(a[counter + 3]))
                    counter +=3
                else:
                    line = j
                o.write(line)
            line_no += 1
if __name__ == "__main__":
    genGro()
    runGromacs()
