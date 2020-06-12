import os, socket

hostname = socket.gethostname()
if 1:
    #os.system('cmdall "tar -xvzf *.tar.gz"')
    os.system('cmdall "cd 2* && python2 ~/soft/simpy/tools/vasp/fe_slow.py "')
    os.system('cat ./*/2*/ave*.dat > fe.dat')
    os.system('python2 ~/soft/simpy/tools/math/int.py')

