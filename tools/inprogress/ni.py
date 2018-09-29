"""
@log:
"""
import numpy as np
import math
import sys
import random

sys.setrecursionlimit(150000)

class Atoms():
    def __init__(self,):
        self.atoms = []
        self.natoms = 0

class Group():
    def __init__(self,):
        self.natoms = 0
        self.atoms = []

def get_occ(dl, r, fname):
    """
    """
    # read the infor
    f = open(fname, "r")

    for i in f: 
        if i.strip().startswith("ITEM: NUMBER OF ATOMS"):
            break
    
    for i in f: 
        if i.strip().startswith("ITEM: BOX BOUNDS"):
            break
        else:
            natom = int(i.strip())
    
    # get the boudary
    boundary = []
    for i in f:
        if i.strip().startswith("ITEM: ATOMS"):
            break
        else:
            boundary.append([float(j) for j in i.strip().split()])
    
    # Grid the space
    xlo = boundary[0][0]
    xhi = boundary[0][1]
    ylo = boundary[1][0]
    yhi = boundary[1][1]
    zlo = boundary[2][0]
    zhi = boundary[2][1]
    nx = int((xhi - xlo)/dl)
    ny = int((yhi - ylo)/dl)
    nz = int((zhi - zlo)/dl)
    dlx = (xhi-xlo)/nx 
    dly = (yhi-ylo)/ny 
    dlz = (zhi-zlo)/nz 
    
    cells = []
    for i in range(nx):
        cells.append([])
        for j in range(ny):
            cells[i].append([])
            for k in range(nz):
                atoms =  Atoms()
                cells[i][j].append(atoms)
    
    for i in f:
        tokens = i.strip().split()
        if len(tokens) > 5:
            n = int(tokens[0])
            id = int(tokens[1])
            x = float(tokens[2])
            y = float(tokens[3])
            z = float(tokens[4])
            charge = 0.0
            x0 = int(round((x - r - xlo)/dlx))
            x1 = int(round((x + r - xlo)/dly))
            y0 = int(round((y - r - ylo)/dly))
            y1 = int(round((y + r - ylo)/dly))
            z0 = int(round((z - r - zlo)/dlz))
            z1 = int(round((z + r - zlo)/dlz))
            for ii in range(x0, x1):
                for jj in range(y0, y1):
                    for kk in range(z0, z1):
                        ii = ii - int(ii/nx)*nx
                        jj = jj - int(jj/ny)*ny
                        kk = kk - int(kk/nz)*nz
                        cells[ii][jj][kk].atoms.append(id)
                        cells[ii][jj][kk].natoms += 1
    f.close()
    
    return cells, nx, ny, nz, dlx, dly, dlz, xlo, ylo, zlo

def cal_occ(cells, nx, ny, nz):
    occ = 0
    un = 0
    reac = 0
    counter = 0

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if cells[i][j][k].natoms > 0:
                    occ += 1
                elif cells[i][j][k].natoms == 0:
                    un += 1
                else:
                    reac += 1
                counter += 1
         
    print "occupied sites = %d with ratio %.4f"%(occ, occ*1.0/counter)
    print "un-occupied sites = ", un*1.0/counter
    print "react sites = %d with ratio %.4f"%( reac, reac*1.0/counter)

def visualize(cells, nx, ny, nz, fname="traj.xyz"):
    """
    """
    n = 0
    o = open(fname, "w")
    o.write("%s\n"%(" "*40))
    o.write("\n")
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if cells[i][j][k].natoms > 0:
                    atp = "O"
                elif cells[i][j][k].natoms == 0:
                    atp = "N"
                else:
                    atp = "C"
                if atp != "N":
                    o.write("%-6s%8d%8d%8d\n"%(atp, i*2, j*2, k*2))
                    n += 1
    o.seek(0)
    o.write("%d"%n)

def get_reac(cells, nx, ny, nz, dlx, dly, dlz, rcut):
    nx_cut = int(rcut/dlx)
    ny_cut = int(rcut/dly)
    nz_cut = int(rcut/dlz)
    print "cut-off for reactions"
    print nx_cut, ny_cut, nz_cut

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if cells[i][j][k].natoms == 0:
                    nv = 0
                    for ii in range(i-nx_cut, i+nx_cut):
                        for jj in range(j-ny_cut, j+ny_cut):
                            for kk in range(k-ny_cut, k+ny_cut):
                                ii = ii - int(ii/nx)*nx
                                jj = jj - int(jj/ny)*ny
                                kk = kk - int(kk/nz)*nz
                                if cells[ii][jj][kk].natoms > 0:
                                    if 3 in cells[ii][jj][kk].atoms:
                                        if 2 in cells[ii][jj][kk].atoms:
                                            nv += 1

                    if nv > 0:
                        cells[i][j][k].natoms = -1

def refine_cells_old(cells, nx, ny):
    """ Delete all the single points
    """

    tmp = []
    for i in range(nx):
        tmp.append([])
        for j in range(ny):
            tmp[i].append(0)

    for i in range(nx):
        for j in range(ny):
            if cells[i][j] == 0:
                nv = 0
                for k in range(i-1, i+2):
                    if k >= nx:
                        k = nx -k
                    for l in range(j-1, j+2):
                        if l >= ny:
                            l = ny - l
                        if cells[k][l] == 0:
                            nv += 1
                if nv == 1:
                    tmp[i][j] = 1
            else:
                tmp[i][j] = cells[i][j]
    return tmp

def refine_cells(cells, grps, ncut):
    """
    """
    for i in grps:
        if i.natoms <= ncut:
            for [j,k,l] in i.atoms:
                cells[j][k][l] = 1


def find_neighbour(cells, nx, ny, nz, x0, y0, z0, visit, boundary):
    nlist  = []

    # x up 
    x = (x0-1) - int((x0-1)/nx)*nx
    y = y0
    z = z0
    nlist.append([x,y,z])
    # x down
    x = (x0+1) - int((x0+1)/nx)*nx
    y = y0
    z = z0
    nlist.append([x,y,z])
    # y up
    x = x0
    y = (y0-1) - int((y0-1)/ny)*ny
    z = z0
    nlist.append([x,y,z])
    # y down
    x = x0
    y = (y0+1) - int((y0+1)/ny)*ny
    z = z0
    # z up
    x = x0
    y = y0
    z = (z0-1) - int((z0-1)/zy)*zy
    nlist.append([x,y,z])
    # z down
    x = x0
    y = y0
    z = (z0+1) - int((z0+1)/nz)*nz
    nlist.append([x,y,z])

    n = 0
    for [i, j, k] in nlist:
        if visit[i][j][k] == 0 and cells[i][j][k] == 1:
            grp.natoms += 1
            visit[i][j][k] = 1
            find_neighbour(cells, nx, ny, nz, i, j, k, visit, grp)
            
    if n > 0 and n < 6:
        boundary.append([nx, ny, nz])
            
def get_clusters(cells, nx, ny, nz):
    grps = []
    visit = []
    for i in range(nx):
        visit.append([])
        for j in range(ny):
            visit[i].append([])
            for k in range(nz):
                visit[i][j].append(0)

    for i in range(nx):
        for j in range(ny):
            if visit[i][j][k] == 0 and cells[i][j][k] == 0:
                grp =  Group()
                grp.natoms = 1
                visit[i][j][k] = 1
                find_neighbour(cells, nx, ny, nz, i, j, k, visit, grp)
                grps.append(grp)
    return grps

def add_reac(cells, nx, ny, dlx, dly, xlo, ylo, fname, atp1, atp2, ratio):
    """
    1 1 0 0 0 0.00901147 0.00607523 -0.000893928
    """
    f = open(fname, "r")
    o = open("replace.lammpstrj", "w")
    for i in f:
        o.write(i)
        if i.startswith("ITEM: ATOMS"):
            break
    for i in f:
        tokens = i.strip().split()
        line = i
        if len(tokens) >5:
            id = int(tokens[1])
            x = float(tokens[2]) - xlo
            y = float(tokens[3]) - ylo
            z = float(tokens[4]) 
            x_n = int(x/dlx)
            y_n = int(y/dly)
            x_n = x_n - int(x_n/nx)*nx
            y_n = y_n - int(y_n/ny)*ny
            if cells[x_n][y_n] == -1 and id == atp1:
                if random.random() < ratio:
                    id = atp2
                tokens[1] = str(id)
                line = " ".join(tokens) + "\n"
        o.write(line)
    o.close()
    f.close()
    
def radius_dist(cells, dlx, dly, grps, ncut):
    """
    """
    ncut = 50
    dist = [0] * 50
    counter = 0
    for i in grps:
        if i.natoms >= ncut:
            counter += 1
            area = i.natoms * dlx * dly
            r = math.sqrt(area/math.pi)
            dist[int(r/2.0)] += 1
    o = open("dist.dat", "w")
    for i in range(len(dist)):
        o.write("%9d%12.4f\n"%(i*2, dist[i]*1.0/counter))
    o.close()
        
    return dist

def add_atoms(cells, nx, ny, nz, dlx, dly, dlz, xlo, ylo, zlo, atp, ratio, infile, outfile):
    lines = []
    f = open(infile, "r")
    for i in f:
        lines.append(i)
    f.close()
    natoms = int(lines[3].strip())

    n_succed = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                if cells[i][j][k].natoms == -1:
                    x = xlo + (i + 0.5) * dlx
                    y = ylo + (j + 0.5) * dly
                    z = zlo + (k + 0.5) * dlz
                    id = atp
                    acc = random.random()
                    if (acc < ratio):
                        n_succed += 1
                        lines.append("%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n"%
                                (natoms + n_succed, id, x, y, z, 0.0, 0.0, 0.0))

    o = open(outfile, "w")
    lines[3] = " %d\n"%(natoms + n_succed)
    for i in lines:
        o.write(i)
    o.close()
    print("Add %d atoms"%n_succed);
            
def main():
    fname = "final.lammpstrj"
    #fname = "test.lammpstrj"
    #dl = 0.936
    dl = 3.0
    r = 1.8 # for Al
    #r = 1.2 # for Ca
    ncut = 2 # cut-off number for cluster size

    cells, nx, ny, nz, dlx, dly, dlz, xlo, ylo, zlo = get_occ(dl, r, fname)
    grps = get_clusters(cells, nx, ny, nz)
    refine_cells(cells, grps, ncut)
    rcut = 4.2 # in A 
    get_reac(cells, nx, ny, nz, dlx, dly, dlz, rcut)
    cal_occ(cells, nx, ny, nz)
    visualize(cells, nx, ny, nz, "reac.xyz")
    atp = 4
    ratio = 0.90
    add_atoms(cells, nx, ny, nz, dlx, dly, dlz, xlo, ylo, zlo, atp, ratio, fname, "add.lammpstrj")

    """
    cal_occ(cells_reac, nx, ny)
    visualize(cells_reac, nx, ny, "reac.xyz")
    atp1 = 3
    atp2 = 4
    ratio = 0.5
    add_reac(cells_reac, nx, ny, dlx, dly, xlo, ylo, fname, atp1, atp2, ratio)
    ncut = 13
    dist = radius_dist(cells, dlx, dly, grps, ncut)
    """

if __name__ == "__main__":
    main()
