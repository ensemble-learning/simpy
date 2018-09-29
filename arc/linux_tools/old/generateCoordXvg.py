def generateCoordXvg(n, N, ndxfile):
    import os
    import os.path
    for i in range(N/n):
        if os.path.exists("m%d"%i):
            pass
        else:
            os.system("mkdir m%d"%i)
        os.system("echo %d | g_traj_mpi -s equil.tpr -f traj.trr -n %s -ox"%((i + 2), ndxfile))
        os.system("mv coord.xvg m%d"%i)
        os.system("cp equil.mdp template.gro *.top *.itp m%d"%i)

if __name__ == "__main__":
    n = 502
    N = 2510
    ndxfile = "out.ndx"
    generateCoordXvg(n, N, ndxfile)    
