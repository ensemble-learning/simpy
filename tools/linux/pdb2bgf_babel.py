import os
import sys

fname = sys.argv[1]

CMD = "babel -ipdb %s.pdb -obgf %s.bgf"%(fname, fname)
os.system(CMD)
