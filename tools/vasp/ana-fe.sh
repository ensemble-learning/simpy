cmdall "tar -xvzf *.tar.gz"
cmdall "cd 4* && python ~/Soft/simpy/tools/vasp/fe_slow.py "
cat ./*/4*/ave*.dat > fe.dat
python ~/Soft/simpy/tools/math/int.py

