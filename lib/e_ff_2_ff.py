import os
from ffield import Ffield
from output_ff import toFfield

ff = Ffield("ffield")
print "Force Field Terms:"
print "atom  : 1"
print "bond  : 2"
print "off   : 3"
print "angle : 4"
term = int(raw_input("Select the term: "))
if term == 1:
    t = ff.atom
    for i in range(len(ff.atom)):
        print i, t[i][0]
elif term == 2:
    t = ff.bond
    for i in range(len(ff.bond)):
        print i, t[i][0], t[i][1]
elif term == 3:
    t = ff.off
    for i in range(len(ff.off)):
        print i, t[i][0], t[i][1]
elif term == 4:
    t = ff.angle
    for i in range(len(ff.angle)):
        print i, t[i][0], t[i][1], t[i][2]
else:
    print "Wrong Choice"
    exit()

nt = int(raw_input("type to change: "))
#temp = float(ff.atom[3][1]) #r0 in eq (2)
print t[nt]

n = int(raw_input("parameter to change: "))

temp = float(ff.bond[nt][n]) #r0 in eq (2)

flag = raw_input("Continue:yes(1);no(0): ")
os.chdir("in")
if int(flag):
    t[nt][n] = "%.4f"%(temp*1.10)
    toFfield(ff, "ffield")
    os.system("bash ./run.sh")
    os.chdir("..")
    os.chdir("de")
    t[nt][n] = "%.4f"%(temp*0.90)
    toFfield(ff, "ffield")
    os.system("./run.sh")
    os.chdir("..")
