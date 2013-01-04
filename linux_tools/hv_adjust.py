# calculate single molecule energy under primary node
# need gromacs input file which can be achieved by dff_2_gmx
# N_TIMES the number of molecule investigated

import os
import os.path
import shutil

LIQUID_DIR = "../"
N_TIMES = 10

def make_input_file(modelFile):
    """ generate gro file use "molecule.gro" as template, and abstract coordinate
        from run.gro. Return the molecule numbers in run.gro, which can be used
        to calculate Hv.
    """
    for nth_times in range(N_TIMES):
        os.system("editconf_mpi	-f %s.gro -o %s.gro -box 5.0 "%(modelFile, modelFile))
        f1 = open("%s.gro"%modelFile, 'r')
        f2 = open("run.gro", 'r')
        o  = open("out_%d.gro"%nth_times, 'w')
        n_atom = 0
        counter = 0
        counter2 = 0
        b = []

        for i in f1:
            if len(i.split()) == 1:
                n_atom = int(i.strip())
                for j in f2:
                    if len(j.split()) == 9 or len(j.split()) == 6:
                        if counter < nth_times*n_atom:
                            pass
                        else:
                            b.append(j[21:].split())
                        counter += 1
                o.write(i)
            elif len(i.split()) == 6:
                a = i[:20]
                o.write("%s%8.3f%8.3f%8.3f\n"%(a,float(b[counter2][0]), float(b[counter2][1]), float(b[counter2][2])))
                counter2 += 1
            else:
                o.write(i)
    return counter / n_atom

def gromacs_job():
    potential_energy = []
    for i in range(N_TIMES):
        os.system("grompp_mpi -f run.mdp -c out_%d.gro -p *.top -o run"%i)
        os.system("mdrun_mpi -s run -e ener_%d.edr"%i)
        read_po_energy = os.popen("echo Potential | g_energy_mpi -f ener_%d.edr"%i)
        for p in read_po_energy:
            if p.strip().startswith("Potential"):
                if p.split()[1].upper() == "NAN":
                    po_energy = 0
                else:
                    po_energy = float(p.split()[1])
                potential_energy.append(po_energy)
        read_po_energy.close()
        for j in os.listdir("./"):
            if j.startswith("#"):
                os.remove(j)
    return potential_energy

def copy_run_gro(target):
    os.system("cp ../*.itp ./")
    potential_energy = 0
    b = os.path.join(LIQUID_DIR, "run.gro")
    ener_file = os.path.join(LIQUID_DIR, "ener.edr")
    if os.path.exists(b):
        shutil.copy(b, target)
        read_po_energy = os.popen("echo Potential | g_energy_mpi -f %s"%ener_file)
        for p in read_po_energy:
            if p.strip().startswith("Potential"):
                potential_energy = float(p.split()[1])
        read_po_energy.close()
        return 1, potential_energy
    else:
        print "file  %s does not exist"%b
        return 0, potential_energy

def single_energy(modelFile):
    potential_energy = {}
    liquid_potential = {}
    molecule_no = {}
    rootDir = os.getcwd()
    copy_sucess, liquid = copy_run_gro(rootDir)
    if copy_sucess:
        molecule_no[modelFile] = make_input_file(modelFile)
        potential_energy[modelFile] = gromacs_job()
        liquid_potential[modelFile] = float(liquid)
    else: 
        pass
    return potential_energy, molecule_no, liquid_potential

def cal_hv(modelFile, temperature):
    #calculation hv
    #This py file should be in the same folder of the single molecule files
    #Two parameters are required: modelFile and temperature
    o = open("hv_result.log", 'w')
    potential_energy, molecule_no, liquid_potential = single_energy(modelFile)
    for i in molecule_no.keys():
        singleEnergy = 0
        for j in potential_energy[i]:
            singleEnergy = singleEnergy + float(j)
        if len(potential_energy[i]) > 0:
            singleEnergy = singleEnergy / len(potential_energy[i])
            hvCal = singleEnergy - liquid_potential[i] / molecule_no[i] + 0.008314 * temperature
            hvCal = hvCal / 4.184
            o.write("%20s%20s%10.3f%4d\n"%(modelFile, i, hvCal, molecule_no[i]))
        else:
            o.write(i + "  error\n")
if __name__ == "__main__":
    cal_hv("dimethyl_ether", 249.9)

