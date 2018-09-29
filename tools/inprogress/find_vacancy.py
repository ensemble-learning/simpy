"""
@log:
"""
import numpy as np
import math
import sys
import random

sys.setrecursionlimit(150000)

def get_occ(dl, r, fname):
    """
    """
    # read the infor
    f = open(fname, "r")

    for i in f: 
        if i.strip().startswith("ITEM: NUMBER OF ATOMS"):
            break
    
    for i in f: 
        if i.strip().startswith("ITEM: BOX BOUNDS pp pp pp"):
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
    nx = int((xhi - xlo)/dl)
    ny = int((yhi - ylo)/dl)
    dlx = (xhi-xlo)/nx 
    dly = (yhi-ylo)/ny 
    
    cells = []
    for i in range(nx):
        cells.append([])
        for j in range(ny):
            cells[i].append(0)
    
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
            for ii in range(x0, x1):
                for jj in range(y0, y1):
                    ii = ii - int(ii/nx)*nx
                    jj = jj - int(jj/ny)*ny
                    cells[ii][jj] += 1
         
    f.close()
    
    return cells, nx, ny, dlx, dly, xlo, ylo

def cal_occ(cells, nx, ny):
    occ = 0
    un = 0
    reac = 0
    counter = 0

    for i in range(nx):
        for j in range(ny):
            if cells[i][j] > 0:
                occ += 1
            elif cells[i][j] == 0:
                un += 1
            else:
                reac += 1
            counter += 1
         
    print "occupied sites = ", occ*1.0/counter
    print "un-occupied sites = ", un*1.0/counter
    print "react sites = ", reac*1.0/counter
    print "react ratio = ", reac*1.0/occ

def visualize(cells, nx, ny, fname="traj.xyz"):
    """
    """
    o = open(fname, "w")
    o.write("%d\n"%(nx * ny))
    o.write("\n")
    for i in range(nx):
        for j in range(ny):
            if cells[i][j] > 0:
                atp = "O"
            elif cells[i][j] == 0:
                atp = "N"
            else:
                atp = "C"
            o.write("%-6s%8d%8d%8d\n"%(atp, i, j, 0))

def get_reac(cells, nx, ny, dlx, dly, rcut):
    nx_cut = int(rcut/dlx)
    ny_cut = int(rcut/dly)
    reac = []
    for i in range(nx):
        reac.append([])
        for j in range(ny):
            reac[i].append(1)
    
    for i in range(nx):
        for j in range(ny):
            if cells[i][j] > 0:
                nv = 0
                for k in range(i-nx_cut, i+nx_cut):
                    if k >= nx:
                        k = nx -k
                    for l in range(j-ny_cut, j+ny_cut):
                        if l >= ny:
                            l = ny - l
                        #r = math.sqrt((k*dlx)**2 + (l*dly)**2)
                        r = 1.0
                        if r < rcut and cells[k][l] == 0:
                            nv += 1
                if nv > 0:
                    reac[i][j] = -1
            else:
                reac[i][j] = 0
    return reac

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
            for [j,k] in i.atoms:
                cells[j][k] = 1

class Group():
    def __init__(self,):
        self.natoms = 0
        self.atoms = []


def find_neighbour(cells, nx, ny, x0, y0, visit, grp):
    nlist  = []

    # up
    x = x0
    y = (y0-1) - int((y0-1)/ny)*ny
    nlist.append([x,y])
    # down
    x = x0
    y = (y0+1) - int((y0+1)/ny)*ny
    nlist.append([x,y])
    # left
    x = (x0-1) - int((x0-1)/nx)*nx
    y = y0
    nlist.append([x,y])
    # right
    x = (x0+1) - int((x0+1)/nx)*nx
    y = y0
    nlist.append([x,y])

    for [i, j] in nlist:
        if visit[i][j] == 0 and cells[i][j] == 0:
            grp.natoms += 1
            grp.atoms.append([i,j])
            visit[i][j] = 1
            find_neighbour(cells, nx, ny, i, j, visit, grp)
        
def get_clusters(cells, nx, ny):
    visit = []
    grps = []
    for i in range(nx):
        visit.append([])
        for j in range(ny):
            visit[i].append(0)

    for i in range(nx):
        for j in range(ny):
            if visit[i][j] == 0 and cells[i][j] == 0:
                grp =  Group()
                grp.natoms = 1
                grp.atoms.append([i,j])
                visit[i][j] = 0
                find_neighbour(cells, nx, ny, i, j, visit, grp)
                grps.append(grp)
    return grps

def test():
    dl = 0.936
    #r = 1.04
    r = 1.2
    cells, nx, ny, dlx, dly = get_occ(dl, r)
    get_clusters(cells, nx, ny)

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
            
def main():
    fname = "final.lammpstrj"
    #fname = "test.lammpstrj"
    dl = 0.936
    r = 1.04 # for Al
    #r = 1.2 # for Ca
    ncut = 10 # cut-off number for cluster size

    cells, nx, ny, dlx, dly, xlo, ylo = get_occ(dl, r, fname)
    grps = get_clusters(cells, nx, ny)
    refine_cells(cells, grps, ncut)
    cal_occ(cells, nx, ny)
    visualize(cells, nx, ny)
    rcut = 12.0 # in A 
    #rcut = 4.5 # in A 
    refine_cells(cells, grps, 5)
    cells_reac = get_reac(cells, nx, ny, dlx, dly, rcut)
    cal_occ(cells_reac, nx, ny)
    visualize(cells_reac, nx, ny, "reac.xyz")
    atp1 = 3
    atp2 = 4
    ratio = 0.25
    add_reac(cells_reac, nx, ny, dlx, dly, xlo, ylo, fname, atp1, atp2, ratio)

    ncut = 13
    dist = radius_dist(cells, dlx, dly, grps, ncut)

if __name__ == "__main__":
    main()
