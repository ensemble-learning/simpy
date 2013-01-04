import os
from ffield import Ffield
from output_ff import toFfield

ff = Ffield("ffield")
temp = float(ff.atom[3][10])
os.chdir("in")
ff.atom[3][10] = "%.4f"%(temp*1.01)
toFfield(ff, "ffield")
os.system("./exe.sh")
os.chdir("..")
os.chdir("de")
ff.atom[3][10] = "%.4f"%(temp*0.99)
toFfield(ff, "ffield")
os.system("./exe.sh")
os.chdir("..")
