import os, shutil

def read_freqs():
    f = open("OUTCAR", "r")

    for i in f:
        if i.strip().startswith("POSCAR ="):
            elements = i.strip().split()[2:]
        if i.strip().startswith("ions per type ="):
            n_atoms = [int(ii) for ii in i.strip().split()[4:]]
        if i.strip().startswith("Eigenvectors and eigenvalues of the dynamical matrix"):
            break

    # begin to read frequency
    freqs = []
    freq_values = []
    flag = 1
    while(flag):
        flag = 0
        for i in f:
            if "cm-1" in i:
                tokens = i.strip().split("=")[1]
                freq_values.append(tokens)
                flag = 1
                freq = []
                break
        for i in f:
            tokens = i.strip().split()
            if len(tokens) == 0:
                freqs.append(freq)
                break
            else:
                freq.append(i)
    f.close()
    return freqs, freq_values

def output_freqs(freqs, freq_values):
    if not os.path.exists("freqs"):
        os.mkdir("freqs")
    os.chdir("freqs")
    for i in range(len(freqs)):
        fn = float(freq_values[i].split()[4])
        folder = "f_%03d_%04d"%(i, fn)
        if not os.path.exists(folder):
            os.mkdir(folder)
        os.chdir(folder)
        o = open("freq", "w")
        o.write("%s\n"%freq_values[i])
        for j in freqs[i]:
            o.write(j)
        o.close()
        os.chdir("..")
    os.chdir("..")

def convert_to_xyz():
    os.chdir("freqs")
    o = open("run.sh", "w")
    o.write("freqmov.pl ../../POSCAR freq 30 0.6")
    o.close()
    os.chdir("..")

freqs, freq_values = read_freqs()
output_freqs(freqs, freq_values)
convert_to_xyz()
