""" Replace the Ca in CaO crystal to creat a Ca-rich LDO (Ca:Al=7:1)
This is for LDH project
"""

from data import ReaxData
from output_conf import toReaxLammps
from tools import LdhProject

for i in range(100):
    testfile = "a666.data"
    a = ReaxData(testfile)
    b = a.parser()
    c = LdhProject(b)
    c.re_al = 8
    c.re1()
    b.assignAtomTypes()
    toReaxLammps(b, "case%02d.data"%i)



