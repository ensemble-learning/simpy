import os
import sys
import shutil
import subprocess
import argparse

def get_2pt(output, tpt="2pt.inp"):
    output.write("2pt analysis configurations:\n")
    
    f = open(tpt, "r")

    for i in f:
        tokens = i.strip().split()
        if "IN_LMPDATA" in i:
            pass
        elif "IN_LMPTRJ" in i:
            pass
        elif "IN_GROUPFILE" in i:
            pass
        elif "ANALYSIS_FRAME_INITIAL" in i:
            pass
        elif "ANALYSIS_FRAME_FINAL" in i:
            pass
        elif "ANALYSIS_FRAME_STEP" in i: 
            pass
        elif "ANALYSIS_VAC_CORLENGTH" in i:
            output.write("    correlation length: %.2f fs\n"%(float(tokens[1])*1000))
        elif "ANALYSIS_VAC_MEMORYMB" in i:
            output.write("    Memory assigned: %d MB\n"%(int(tokens[1])))
        elif "ANALYSIS_VAC_2PT" in i:
            pass
        elif "ANALYSIS_VAC_FIXED_DF" in i:
            output.write("    Remove %d freedoms\n"%(int(tokens[1])))
        elif "ANALYSIS_OUT" in i:
            pass
        elif "ANALYSIS_LMP_TSTEP" in i:
            output.write("    Lammps timestep %.2f fs\n"%(float(tokens[1])*1000))
            pass
        elif "ANALYSIS_VAC_LINEAR_MOL" in i:
            pass
        elif "ANALYSIS_VAC_ROTN_SYMMETRY" in i:
            pass

def get_lammps_inp(output, lammps_input="lammps_input"):

    unit = 0
    timestep = 0.0
    savefreq = 0
    steps = 0

    f = open(lammps_input, "r")
    for i in f:
        tokens = i.strip().split()
        if i.strip().startswith("#"):
            pass
        elif len(tokens) == 0:
            pass
        else:
            if "units" in i:
                if tokens == "real":
                    unit = 1
            elif "timestep" in i:
                timestep = float(tokens[1])
            elif "dump" == tokens[0]:
                savefreq = float(tokens[4])
            elif "run" == tokens[0]:
                steps = int(tokens[1]) 
    output.write("Lammps infor\n")
    output.write("    timestep = %6.2f fs\n"%(timestep))
    output.write("    total = %6.2f ps\n"%(steps * timestep /1000.0))
    output.write("    save trj every %6.2f fs\n"%(timestep * savefreq))

def run_2pt():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", nargs=1, type=int, help="number of CPU")
    parser.add_argument("-t", nargs=1, type=int, help="number of 2pt iterations")
    args = parser.parse_args()

    o = open("run.log", "w")

    get_lammps_inp(o)
    o.flush()
    
    get_2pt(o)

    if args.n:
        ppn = args.n[0]
    else:
        ppn = 16

    if args.t:
        iters = args.t[0]
    else:
        iters = 5
    mpirun = "/net/hulk/home6/chengtao/bin/openmpi/bin/mpirun"
    
    o.write("Using %d cpu and %d steps\n"%(ppn, iters))
    o.flush()
    
    for i in range(iters):
        o.write("Runing step %8d\n"%i)
        # prepare input files
        folder = "sample_%02d"%i
        if not os.path.exists(folder):
            os.mkdir(folder)
        shutil.copy("lammps_input", folder)
        shutil.copy("lammps_restart", folder)
        shutil.copy("lammps.data", folder)
        shutil.copy("2pt.inp", folder)
        # do the calculation
        os.chdir(folder)
        log = open("lammps.run", "w")
        error = open("lammps.err", "w")
        o.write("    Runing Lammps\n")
        signal = subprocess.call([mpirun, "-np", str(ppn), "lammps", "-in", "lammps_input"], \
                 stdout=log, stderr=error)
        if signal > 0:
            sys.exit()
        o.write("    Normal end Lammps\n")
        log.close()
        error.close()
        shutil.copy("log.lammps", "inp.eng")
        shutil.move("dump.lammpstrj", "inp.lammps")
        log = open("2pt.run", "w")
        error = open("2pt.err", "w")
        o.write("    Runing 2pt\n")
        subprocess.call(["2pt_analysis", "2pt.inp"], stdout=log, stderr=error)
        if signal > 0:
            sys.exit()
        o.write("    Normal end 2pt\n")
        log.close()
        error.close()
        os.remove("inp.lammps")
        os.chdir("..")
        shutil.copy(os.path.join(folder,"lammps_conti"), "lammps_restart")
        o.flush()
    o.close()

if __name__ == "__main__":
    run_2pt()
