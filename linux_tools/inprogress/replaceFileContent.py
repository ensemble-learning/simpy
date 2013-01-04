import os
import sys
import string

def gromacs():
        os.system('grompp_mpi -np 2 -f equil -c confout.gro -p *.top -o equil')
        os.system('mpirun -np 2 mdrun_mpi -s equil')

def replace( file, search_for, replace_with):
    # replace strings in a text file

    back = file + ".bak"
    temp = file + ".tmp"

    try:
        # remove old temp file, if any
        os.remove(temp)
    except os.error:
        pass

    fi = open(file)
    fo = open(temp, "w")

    for s in fi.readlines():
        fo.write(string.replace(s, search_for, replace_with))

    fi.close()
    fo.close()

    try:
        # remove old backup file, if any
        os.remove(back)
    except os.error:
        pass

    # rename original to backup...
    os.rename(file, back)

    # ...and temporary to original
    os.rename(temp, file)

#
# try it out!

if __name__ == '__main__':
    iteration = 3
    for i in range(iteration):
        for j in range(273,423,25):
            os.chdir('./%d'%j)
            replace('towhee_input', '/home/hsun/xfli/oxide/towhee_ff_oxide4', './towhee_ff_oxide4')
            #gromacs()
            os.chdir(os.pardir)

