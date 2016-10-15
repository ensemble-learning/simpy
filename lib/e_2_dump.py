from block import dumpBlock
from dump import Dump
from output_conf import toXyz, toPdb, toPoscar

lmpfile = "dump.lammpstrj"
sepfile = "dump.sep"
dt = 1


# parse the dump file with multi configurations into seperated dump files
nframe = dumpBlock(lmpfile, sepfile, dt)
nframe += 1

for i in range(0, nframe, dt):
#for i in range(10):
    a = Dump("%s%05d.dump"%(sepfile,i))
    b = a.parser()
    b.assignEleTypes()
    b.toFrac()

    toXyz(b, "xyz%05d.xyz"%i)
    #b.sortXYZ("z")
    toPdb(b, "pdb%05d.pdb"%i)
    toReaxLammps(b, "lammps.data")
    toPoscar(b,)




