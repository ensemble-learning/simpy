import os
from g03 import G03tools

os.chdir("d:\\Nut1\\work\\cement\\crystal_new\\Belite\\reactions\\scan\\a9\\")
for i in range(10):
    b = G03tools('scan%d.log'%i)
    print b.getEnergy()



