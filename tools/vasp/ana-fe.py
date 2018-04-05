import os, socket

hostname = socket.gethostname()
if "stampede2.tacc.utexas.edu" in hostname:
    os.system('cmdall "tar -xvzf *.tar.gz"')
    os.system('cmdall "cd 6* && python ~/soft/simpy/tools/vasp/fe_slow.py "')
    os.system('cat ./*/6*/ave*.dat > fe.dat')
    os.system('python ~/soft/simpy/tools/math/int.py')

