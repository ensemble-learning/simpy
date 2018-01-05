cmdall "tar -xvzf *.tar.gz"
cmdall "cd 5* && python ~/Soft/simpy/tools/vasp/fe_slow.py "
cat ./*/5*/ave*.dat > fe.dat
python ~/Soft/simpy/tools/math/int.py

