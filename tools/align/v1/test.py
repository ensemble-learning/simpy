import math
import numpy as np
import scipy.linalg

def read_coords():
    f = open("test.xyz", "r")
    lines = f.readlines()
    coords = []
    for i in range(2, len(lines)):
        tokens = lines[i].strip().split()
        if len(tokens) > 0:
            coords.append([float(j) for j in tokens[1:]])
    return coords

def get_surface_atoms(coords):
    sur = []
    for i in range(32):
        sur.append(coords[i])
    return sur

def get_com(coords):
    com = [0.0, 0.0, 0.0]
    for i in range(len(coords)):
        for j in range(3):
            com[j] += coords[i][j]
    for i in range(3):
        com[i] = com[i]/len(coords)

    return com

def center_coords(coords, s):
    coords_new = []
    for i in range(len(coords)):
        coords_new.append([])
        for j in range(3):
            coords_new[i].append(coords[i][j] - s[j])
    return coords_new
        
def fit_plannar(coords):
    # regular grid covering the domain of the data
    debug = 0
    data = np.array(coords)
    X,Y = np.meshgrid(np.arange(-10.0, 10.0, 1), np.arange(-10.0, 10.0, 1))
    # best-fit linear plane
    A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients
    if debug:
        o = open("plannar.xyz", "w")
        o.write("%d\n\n"%(len(coords)*2))
        for i in range(len(coords)):
            x = coords[i][0]
            y = coords[i][1]
            z0 = coords[i][2]
            z1 = C[0]*x + C[1]*y + C[2]
            o.write("Au    %12.4f%12.4f%12.4f\n"%(x, y, z0))
            o.write("Fe    %12.4f%12.4f%12.4f\n"%(x, y, z1))
        o.close()
    return C

def align_axis(v1, v2):
    w = np.cross(v1,v2)
    w = w/np.linalg.norm(w)
    w_hat = np.matrix([[    0, -w[2],  w[1]],
                                              [ w[2],     0, -w[0]],
                                              [-w[1],  w[0],     0]])
    cos_tht = np.dot(v1, v2)/np.linalg.norm(v1)/np.linalg.norm(v2)
    tht = np.arccos(cos_tht)
    R = np.identity(3) + w_hat*np.sin(tht) + w_hat*w_hat*(1-np.cos(tht))
    return R

def rotate_coords(coords, R):
    coords_new = []
    for i in range(len(coords)):
        coords_new.append(np.dot(coords[i], R).A1)
    return coords_new

def output_coords(coords, fname):
    o = open(fname, "w")
    o.write("%d\n\n"%(len(coords)+20))
    for i in range(len(coords)):
        o.write("Au    ")
        for j in range(3):
            o.write("%12.4f"%coords[i][j])
        o.write("\n")
    for i in range(20):
        o.write("O 0 0 %12.4f\n"%(i-10))
    o.close()

def output_coords_debug(coords, fname, pc):
    o = open(fname, "w")
    o.write("%d\n\n"%(len(coords)+20))
    for i in range(len(coords)):
        o.write("Au    ")
        for j in range(3):
            o.write("%12.4f"%coords[i][j])
        o.write("\n")
    for i in range(20):
        x = pc[0]*i
        y = pc[1]*i
        z = -i
        o.write("O %12.4f%12.4f%12.4f\n"%(x, y, z))
    o.close()

debug = 0
# read the xyz file
coords = read_coords()
# get the surface atoms
sur = get_surface_atoms(coords)
# calculate the center of mass (com)
com = get_com(sur)
# move the com to zero
sur_zero = center_coords(sur, com)
# determine the surface and surface vector
pc = fit_plannar(sur_zero)
v1 = [pc[0], pc[1], -1]
v2 = [0, 0, 1]
# calculate the rotation matrix
R = align_axis(v1, v2)
print np.dot(R, v1)
# apply the rotation
coords_zero = center_coords(coords, com)
coords_new = rotate_coords(coords_zero, R)
# output the data
output_coords(coords_zero, "t0.xyz")
output_coords(coords_zero, "t0.xyz")
if debug:
    output_coords_debug(coords_zero, "vector.xyz", pc)
