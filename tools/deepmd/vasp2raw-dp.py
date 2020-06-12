from dpdata import System, LabeledSystem, MultiSystems

#s = System('POSCAR', fmt='poscar')
#print(s)
ls=LabeledSystem('OUTCAR', fmt='outcar')
"""
if len(ls)%2==0:
    size = int(len(ls)/2)
else:
    size = int(len(ls)/2) + 1
"""
ls.to_deepmd_raw('.')
#ls.to_deepmd_npy('deepmd', set_size=size)
