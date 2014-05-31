import os

FOLDERS = ["BCC", "DIAMOND", "FCC", "HCP", "SCC"]
ELEMENT = "Ti"
FLAGS = []
for i in FOLDERS:
    FLAGS.append(ELEMENT + i)

get_energy = "bash ~/soft/simpy/tools/vasp/sum_energy.sh"
get_vols = "cmdall 'python ~/soft/simpy/lib/e_2_contcar.py' > vols"
get_bulk_modulus = "python ~/soft/simpy/tools/vasp/eos.py"
get_geo = "cat ./*/geo > geo"
add_flag = "python ~/soft/simpy/tools/reax/addflag.py %flag%"

n = 0
for i in FOLDERS:
    os.chdir(i)
    os.chdir("scan")
    os.system(get_energy)
    os.system(get_vols)
    #os.system(get_bulk_modulus)
    #os.system(get_geo)
    #cmd = add_flag
    #cmd = cmd.replace("%flag%", FLAGS[n])
    #os.system(cmd)
    os.chdir("..")
    os.chdir("..")
    n += 1

