from block import dumpBlock
from dump import Dump
from output_conf import toXyz, toPdb

# parse the dump file with multi configurations into seperated dump files
dumpBlock("dump.lmp")

for i in range(200):
    a = Dump("dump.sep%05d.dump"%i)
    b = a.parser()
    toXyz(b, "xyz%02d.xyz"%i)
    toPdb(b, "pdb%02d.pdb"%i)




