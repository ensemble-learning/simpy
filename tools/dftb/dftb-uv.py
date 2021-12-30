import sys, os
import shutil

dpath = '/home/tcheng/soft/simpy/tools/dftb'

fname = sys.argv[1]
os.system('python %s %s'%(os.path.join(dpath, 'any2gen'), fname))
os.system('python %s %s'%(os.path.join(dpath, 'get_input.py'), 'opt'))
os.system('dftb+')
shutil.copy('geo_end.gen', 'geo_start.gen')
os.system('python %s %s'%(os.path.join(dpath, 'get_input.py'), 'tddft'))
os.system('dftb+')
os.system('python %s'%(os.path.join(dpath, 'uv.py')))
