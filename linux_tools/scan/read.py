import os
os.system("ls ./*/*.edr > ls.log")
os.system("echo '15' > in.log")
os.system("echo '11' >> in.log")
os.system("echo '14' >> in.log")
os.system("echo '20' >> in.log")
f = open("ls.log")
for i in f:
    a = i[:-1]
    os.system("echo %s >> out.log"%a)
    os.system("g_energy_mpi -f %s < in.log >> out.log"%a)
#os.system("grep './' > list.log")
#os.system("grep 'dVpot/dlambda' > result.log")

