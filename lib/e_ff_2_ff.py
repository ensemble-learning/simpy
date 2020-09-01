import os
from ffield import Ffield
from output_ff import toFfield

ff = Ffield("ffield", 0)
print("Force Field Terms:")
print("atom  : 1")
print("bond  : 2")
print("off   : 3")
print("angle : 4")
term = int(input("Select the term: "))
if term == 1:
    t = ff.atom
    for i in range(len(ff.atom)):
        print(i, t[i][0])
elif term == 2:
    t = ff.bond
    for i in range(len(ff.bond)):
        print(i, t[i][0], t[i][1])
elif term == 3:
    t = ff.off
    for i in range(len(ff.off)):
        print(i, t[i][0], t[i][1])
elif term == 4:
    t = ff.angle
    for i in range(len(ff.angle)):
        print(i, t[i][0], t[i][1], t[i][2])
else:
    print("Wrong Choice")
    exit()

nt = int(input("type to change: "))
#temp = float(ff.atom[3][1]) #r0 in eq (2)
print(t[nt])

n = int(input("parameter to change: "))

temp = float(t[nt][n]) #r0 in eq (2)

scale = float(input("scale factor: "))

flag = input("Continue:yes(1);no(0): ")
if int(flag):
    # Increase
    if not os.path.exists("in"):
        os.mkdir("in")
    os.chdir("in")
    increase = 1.0 + scale
    t[nt][n] = "%.4f"%(temp*increase)
    toFfield(ff, "ffield")
    os.system("bash ./run.sh")
    os.chdir("..")
    # Original
    if not os.path.exists("or"):
        os.mkdir("or")
    os.chdir("or")
    origin = 1.0
    t[nt][n] = "%.4f"%(temp*origin)
    toFfield(ff, "ffield")
    os.system("bash ./run.sh")
    os.chdir("..")
    # Decrease
    if not os.path.exists("de"):
        os.mkdir("de")
    os.chdir("de")
    decrease = 1.0 - scale
    t[nt][n] = "%.4f"%(temp*decrease)
    toFfield(ff, "ffield")
    os.system("bash ./run.sh")
    os.chdir("..")
