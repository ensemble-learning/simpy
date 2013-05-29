# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under 
# the GNU General Public License.

# patch tool

oneline = "Create patchy Lennard-Jones particles for LAMMPS input"

docstr = """
p = patch(vfrac)           setup box with a specified volume fraction
p = patch(vfrac,1,1,2)     x,y,z = aspect ratio of box (def = 1,1,1)
    
p.seed = 48379		   set random # seed (def = 12345)
p.dim = 2		   set dimension of created box (def = 3)
p.blen = 0.97              set length of tether bonds (def = 0.97)
p.dmin = 1.02              set min r from i-1 to i+1 tether site (def = 1.02)

p.build(100,"hex2",1,2,3)  create 100 "hex2" particles with types 1,2,3
  
  can be invoked multiple times
  keywords:
    c60hex2: diam,1,2,3 = C-60 with 2 hex patches and ctr part, types 1,2,3
    hex2: diam,1,2 = one large particle with 2 7-mer hex patches, types 1,2
    hex4: diam,1,2 = one large particle with 4 7-mer hex patches, types 1,2
    ring: diam,N,1,2 = one large part with equatorial ring of N, types 1,2
    ball: diam,m1,m2,1,2,3 = large ball with m12-len tethers, types 1,2,3
    tri5: 1,2 = 3-layer 5-size hollow tri, types 1,2
    rod: N,m1,m2,1,2,3 = N-length rod with m12-len tethers, types 1,2,3
    tri: N,m1,m2,m3,1,2,3,4 = N-size tri with m123-len tethers, types 1-4
    hex: m1,m2,m3,m4,m5,m6,1,2,3,4,5,6,7 = 7-atom hex with m-len tethers, t 1-7
    dimer: r,1 = two particles r apart, type 1, no bond
    star2d: N,r,1 = 2d star of length N (odd), beads r apart, type 1, no bonds
    box2d: N,M,r,1 = 2d NxM box, beads r apart, type 1, no bonds
    
p.write("data.patch")      write out system to LAMMPS data file
"""

# History
#   8/05, Steve Plimpton (SNL): original version

# ToDo list

# Variables
#   vfrac = desired volume fraction
#   x,y,z = aspect ratio of box (def = 1,1,1)
#   seed = random seed
#   molecules = list of atoms 

# Imports and external programs

from math import sqrt,pi,cos,sin
from data import data

# Class definition

class patch:
  
  # --------------------------------------------------------------------

  def __init__(self,vfrac,*list):
    self.vfrac = vfrac
    self.xaspect = self.yaspect = self.zaspect = 1.0
    if len(list):
      self.xaspect = list[0]
      self.yaspect = list[1]
      self.zaspect = list[2]
    self.seed = 12345
    self.dim = 3
    self.molecules = []
    self.volume = 0
    self.blen = 0.97
    self.dmin = 1.02

  # --------------------------------------------------------------------
  # call style method with extra args
  # adds to volume and atom list

  def build(self,n,style,*types):
    cmd = "atoms,bonds,volume = self.%s(*types)" % style
    for i in xrange(n):
      exec cmd
      self.molecules.append([atoms,bonds])
      self.volume += volume

  # --------------------------------------------------------------------
  # create the atom coords in a scaled-size box of correct dimension
  # write them to LAMMPS data file

  def write(self,file):
    if self.dim == 3: self.write3d(file)
    else: self.write2d(file)

  # --------------------------------------------------------------------
  # write a 3d simulation to data file

  def write3d(self,file):
    volume = self.volume/self.vfrac
    prd = pow(volume/self.xaspect/self.yaspect/self.zaspect,1.0/3.0)
    self.xprd = self.xaspect * prd
    self.xlo = -self.xprd/2.0
    self.xhi = self.xprd/2.0
    self.yprd = self.yaspect * prd
    self.ylo = -self.yprd/2.0
    self.yhi = self.yprd/2.0
    self.zprd = self.zaspect * prd
    self.zlo = -self.zprd/2.0
    self.zhi = self.zprd/2.0

    idatom = idbond = idmol = 0
    atoms = []
    bonds = []
    xp = 3*[0]
    yp = 3*[0]
    zp = 3*[0]

    for molecule in self.molecules:
      idmol += 1

      # xp[3],yp[3],zp[3] = randomly oriented, normalized basis vectors
      # xp is in random direction
      # yp is random dir crossed into xp
      # zp is xp crossed into yp

      xp[0] = self.random() - 0.5
      xp[1] = self.random() - 0.5
      xp[2] = self.random() - 0.5
      r = sqrt(xp[0]*xp[0] + xp[1]*xp[1] + xp[2]*xp[2])
      xp[0],xp[1],xp[2] = xp[0]/r,xp[1]/r,xp[2]/r

      r0 = self.random() - 0.5
      r1 = self.random() - 0.5
      r2 = self.random() - 0.5
      yp[0] = r1*xp[2] - r2*xp[1]
      yp[1] = r2*xp[0] - r0*xp[2]
      yp[2] = r0*xp[1] - r1*xp[0]
      r = sqrt(yp[0]*yp[0] + yp[1]*yp[1] + yp[2]*yp[2])
      yp[0],yp[1],yp[2] = yp[0]/r,yp[1]/r,yp[2]/r

      zp[0] = xp[1]*yp[2] - xp[2]*yp[1]
      zp[1] = xp[2]*yp[0] - xp[0]*yp[2]
      zp[2] = xp[0]*yp[1] - xp[1]*yp[0]
      r = sqrt(zp[0]*zp[0] + zp[1]*zp[1] + zp[2]*zp[2])
      zp[0],zp[1],zp[2] = zp[0]/r,zp[1]/r,zp[2]/r

      # random origin for new particle

      xorig = self.xlo + self.random()*self.xprd
      yorig = self.ylo + self.random()*self.yprd
      zorig = self.zlo + self.random()*self.zprd

      # unpack bonds in molecule before atoms so idatom = all previous atoms
      
      for bond in molecule[1]:
        idbond += 1
        bonds.append([idbond,bond[0],bond[1]+idatom+1,bond[2]+idatom+1])

      # unpack atoms in molecule
      # xnew,ynew,xnew = coeffs in new rotated basis vectors
      # x,y,z = coeffs in original xyz axes

      for atom in molecule[0]:
        idatom += 1
	xnew = atom[1]
	ynew = atom[2]
	znew = atom[3]
        x = xorig + xnew*xp[0] + ynew*yp[0] + znew*zp[0]
        y = yorig + xnew*xp[1] + ynew*yp[1] + znew*zp[1]
        z = zorig + xnew*xp[2] + ynew*yp[2] + znew*zp[2]
	ix = iy = iz = 0
	x,y,z,ix,iy,iz = self.pbc(x,y,z,ix,iy,iz)
        atoms.append([idatom,idmol,atom[0],x,y,z,ix,iy,iz])

    list = [atom[2] for atom in atoms]
    atypes = max(list)

    # create the data file

    d = data()
    d.title = "LAMMPS data file for Nanoparticles"
    d.headers["atoms"] = len(atoms)
    d.headers["atom types"] = atypes
    if bonds:
      d.headers["bonds"] = len(bonds)
      d.headers["bond types"] = 1
    d.headers["xlo xhi"] = (self.xlo,self.xhi)
    d.headers["ylo yhi"] = (self.ylo,self.yhi)
    d.headers["zlo zhi"] = (self.zlo,self.zhi)

    lines = []
    for atom in atoms:
      line = "%d %d %d %g %g %g %d %d %d\n" % \
             (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5],
              atom[6], atom[7], atom[8])
      lines.append(line)
    d.sections["Atoms"] = lines

    if bonds:
      lines = []
      for bond in bonds:
        line = "%d %d %d %d\n" % (bond[0], bond[1], bond[2], bond[3])
        lines.append(line)
        d.sections["Bonds"] = lines
      
    d.write(file)

  # --------------------------------------------------------------------
  # write a 2d simulation to data file

  def write2d(self,file):
    volume = self.volume/self.vfrac
    prd = pow(volume/self.xaspect/self.yaspect,1.0/2.0)
    self.xprd = self.xaspect * prd
    self.xlo = -self.xprd/2.0
    self.xhi = self.xprd/2.0
    self.yprd = self.yaspect * prd
    self.ylo = -self.yprd/2.0
    self.yhi = self.yprd/2.0
    self.zlo = -0.5
    self.zhi = 0.5

    idatom = idbond = idmol = 0
    atoms = []
    bonds = []
    xp = 3*[0]
    yp = 3*[0]

    for molecule in self.molecules:
      idmol += 1

      # xp[2],yp[2] = randomly oriented, normalized basis vectors
      # xp is in random direction
      # yp is (0,0,1) crossed into xp

      xp[0] = self.random() - 0.5
      xp[1] = self.random() - 0.5
      r = sqrt(xp[0]*xp[0] + xp[1]*xp[1])
      xp[0],xp[1] = xp[0]/r,xp[1]/r

      yp[0] = -xp[1]
      yp[1] = xp[0]
      r = sqrt(yp[0]*yp[0] + yp[1]*yp[1])
      yp[0],yp[1] = yp[0]/r,yp[1]/r

      # random origin for new particle

      xorig = self.xlo + self.random()*self.xprd
      yorig = self.ylo + self.random()*self.yprd
      zorig = 0.0

      # unpack bonds in molecule before atoms so idatom = all previous atoms
      
      for bond in molecule[1]:
        idbond += 1
        bonds.append([idbond,bond[0],bond[1]+idatom+1,bond[2]+idatom+1])

      # unpack atoms in molecule
      # xnew,ynew,xnew = coeffs in new rotated basis vectors
      # x,y,z = coeffs in original xyz axes

      for atom in molecule[0]:
        idatom += 1
	xnew = atom[1]
	ynew = atom[2]
        x = xorig + xnew*xp[0] + ynew*yp[0]
        y = yorig + xnew*xp[1] + ynew*yp[1]
        z = 0.0
	ix = iy = iz = 0
	x,y,z,ix,iy,iz = self.pbc(x,y,z,ix,iy,iz)
        atoms.append([idatom,idmol,atom[0],x,y,z,ix,iy,iz])

    list = [atom[2] for atom in atoms]
    atypes = max(list)

    # create the data file

    d = data()
    d.title = "LAMMPS data file for Nanoparticles"
    d.headers["atoms"] = len(atoms)
    d.headers["atom types"] = atypes
    if bonds:
      d.headers["bonds"] = len(bonds)
      d.headers["bond types"] = 1
    d.headers["xlo xhi"] = (self.xlo,self.xhi)
    d.headers["ylo yhi"] = (self.ylo,self.yhi)
    d.headers["zlo zhi"] = (self.zlo,self.zhi)

    lines = []
    for atom in atoms:
      line = "%d %d %d %g %g %g %d %d %d\n" % \
             (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5],
              atom[6], atom[7], atom[8])
      lines.append(line)
    d.sections["Atoms"] = lines

    if bonds:
      lines = []
      for bond in bonds:
        line = "%d %d %d %d\n" % (bond[0], bond[1], bond[2], bond[3])
        lines.append(line)
        d.sections["Bonds"] = lines
      
    d.write(file)

  # --------------------------------------------------------------------

  def pbc(self,x,y,z,ix,iy,iz):
    if x < self.xlo:
      x += self.xprd
      ix -= 1
    elif x >= self.xhi:
      x -= self.xprd
      ix += 1
    if y < self.ylo:
      y += self.yprd
      iy -= 1
    elif y >= self.yhi:
      y -= self.yprd
      iy += 1
    if z < self.zlo:
      z += self.zprd
      iz -= 1
    elif z >= self.zhi:
      z -= self.zprd
      iz += 1
    return x,y,z,ix,iy,iz

  # --------------------------------------------------------------------
  # params = diam,type1,type2,type3
  # type1 = type of non-patch atoms, type2 = type of patch atoms
  # type3 = type of center-of-sphere atom
  
  def c60hex2(self,*params):
    template = BUCKY_60
    diam = params[0]
    patches = [[params[2],[1,2,3,4,5,6,34,35,46,47,48,60]],
               [params[3],[61]]]
    atoms = make_sphere(template,diam,params[1],patches)
    volume = 4.0/3.0 * pi * diam*diam*diam/8
    return atoms,[],volume
  
  # --------------------------------------------------------------------
  # params = diam,type1,type2
  # type1 = type of large center atom, type2 = type of hex patch atoms
  
  def hex2(self,*params):
    diam = params[0]
    type1 = params[1]
    type2 = params[2]
    
    atoms = []
    atoms.append([type1,0.0,0.0,0.0])
    
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,0.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,-1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,-0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,0.5,-sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,-0.5,-sqrt(3.0)/2))

    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,0.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,-1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,-0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,0.5,-sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,-0.5,-sqrt(3.0)/2))

    volume = 4.0/3.0 * pi * diam*diam*diam/8
    return atoms,[],volume
  
  # --------------------------------------------------------------------
  # params = diam,type1,type2
  # type1 = type of large center atom, type2 = type of hex patch atoms
  
  def hex4(self,*params):
    diam = params[0]
    type1 = params[1]
    type2 = params[2]
    
    atoms = []
    atoms.append([type1,0.0,0.0,0.0])
    
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,0.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,-1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,-0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,0.5,-sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5*diam,-0.5,-sqrt(3.0)/2))

    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,0.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,-1.0,0.0))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,-0.5,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,0.5,-sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5*diam,-0.5,-sqrt(3.0)/2))

    atoms.append(atom_on_sphere(diam,type2,0.0,0.5*diam,0.0))
    atoms.append(atom_on_sphere(diam,type2,1.0,0.5*diam,0.0))
    atoms.append(atom_on_sphere(diam,type2,-1.0,0.5*diam,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5,0.5*diam,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5,0.5*diam,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5,0.5*diam,-sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5,0.5*diam,-sqrt(3.0)/2))

    atoms.append(atom_on_sphere(diam,type2,0.0,-0.5*diam,0.0))
    atoms.append(atom_on_sphere(diam,type2,1.0,-0.5*diam,0.0))
    atoms.append(atom_on_sphere(diam,type2,-1.0,-0.5*diam,0.0))
    atoms.append(atom_on_sphere(diam,type2,0.5,-0.5*diam,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5,-0.5*diam,sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,0.5,-0.5*diam,-sqrt(3.0)/2))
    atoms.append(atom_on_sphere(diam,type2,-0.5,-0.5*diam,-sqrt(3.0)/2))

    volume = 4.0/3.0 * pi * diam*diam*diam/8
    return atoms,[],volume

  # --------------------------------------------------------------------
  # params = diam,nring,type1,type2
  # nring = # of particles in ring
  # type1 = type of large center atom, type2 = type of ring atoms
  
  def ring(self,*params):
    diam = params[0]
    nring = params[1]
    type1 = params[2]
    type2 = params[3]
    
    atoms = []
    atoms.append([type1,0.0,0.0,0.0])

    for i in range(nring):
      atoms.append([type2,0.5*diam*cos(i * 2*pi/nring),
                    0.5*diam*sin(i * 2*pi/nring),0.0])

    volume = 4.0/3.0 * pi * diam*diam*diam/8
    return atoms,[],volume

  # --------------------------------------------------------------------
  # params = diam,m1,m2,ntype,m1type,m2type
  # ntype = type of big central ball
  # m12,m12type = length of tethers on each side of ball (m12 = 0 = no tether)
  # set three types of bonds:
  #   1 = big to small, 2 = small to small, 3 = across two tethers
  
  def ball(self,*params):
    diam = params[0]
    m1 = params[1]
    m2 = params[2]
    ntype = params[3]
    m1type = params[4]
    m2type = params[5]

    atoms = []
    atoms.append([ntype,0.0,0.0,0.0])

    if m1:
      atoms.append([m1type,0.5*diam+0.5,0.0,0.0])
      atoms += tether(m1-1,m1type,self.blen,self.dmin,
                      [ntype,0.5*diam-0.5,0.0,0.0],atoms[1],self.random)
    if m2:
      atoms.append([m2type,-0.5*diam-0.5,0.0,0.0])
      atoms += tether(m2-1,m2type,self.blen,self.dmin,
                      [ntype,-0.5*diam+0.5,0.0,0.0],atoms[1+m1],self.random)

    bonds = []
    for i in range(m1):
      if i == 0: bonds.append([1,0,1])
      else: bonds.append([2,1+i-1,1+i])
    for i in range(m2):
      if i == 0: bonds.append([1,0,1+m1])
      else: bonds.append([2,1+m1+i-1,1+m1+i])
    if m1 and m2: bonds.append([3,1,m1+1])
      
    volume = 4.0/3.0 * pi * diam*diam*diam/8 + (m1+m2)*pi/6.0
    return atoms,bonds,volume

  # --------------------------------------------------------------------
  # params = type1,type2
  # type1 = type of 1/3 layers, type2 = type of middle layer
  
  def tri5(self,*params):
    template = TRI5_HOLLOW
    nlayer = 3
    patches = [[params[1],[13,14,15,16,17,18,19,20,21,22,23,24]]]

    atoms = []
    n = 0
    for iz in range(nlayer):
      for atom in template:
        n += 1
        type = params[0]
        for pair in patches:
          patch = pair[1]
          if n in patch: type = pair[0]
        x,y,z = atom[0],atom[1],atom[2] + iz
        atoms.append([type,x,y,z])

    volume = nlayer * 0.5 * 5.0 * 5.0*sqrt(3)/2
    return atoms,[],volume

  # --------------------------------------------------------------------
  # params = n,m1,m2,ntype,m1type,m2type
  # n,ntype = length/type of rod
  # m12,m12type = length of tethers on each end (m12 = 0 = no tether)

  def rod(self,*params):
    n = params[0]
    m1 = params[1]
    m2 = params[2]
    ntype = params[3]
    m1type = params[4]
    m2type = params[5]
    
    atoms = []
    for i in range(n):
      x,y,z = i*self.blen,0.0,0.0
      atoms.append([ntype,x,y,z])
      
    if m1: atoms += tether(m1,m1type,self.blen,self.dmin,
                           atoms[n-2],atoms[n-1],self.random)
    if m2: atoms += tether(m2,m2type,self.blen,self.dmin,
                           atoms[1],atoms[0],self.random)

    bonds = []
    for i in range(m1):
      if i == 0: bonds.append([1,n-1,n])
      else: bonds.append([1,n+i-1,n+i])
    for i in range(m2):
      if i == 0: bonds.append([1,0,n+m1])
      else: bonds.append([1,n+m1+i-1,n+m1+i])
  
    volume = (n+m1+m2) * pi / 6.0
    return atoms,bonds,volume

  # --------------------------------------------------------------------
  # params = nsize,m1,m2,m3,ntype,m1type,m2type,m3type
  # nsize,ntype = size,type of triangle center
  # m123,m123type = length of tethers on each corner (m123 = 0 = no tether)
  
  def tri(self,*params):
    nsize = params[0]
    m1 = params[1]
    m2 = params[2]
    m3 = params[3]
    ntype = params[4]
    m1type = params[5]
    m2type = params[6]
    m3type = params[7]
  
    atoms = []
    for i in range(nsize):
      n = nsize - i
      for j in range(n):
        x,y,z = 0.5*i*self.blen + j*self.blen, i*self.blen*sqrt(3.0)/2.0, 0.0
        atoms.append([ntype,x,y,z])

    n = len(atoms)
    if m1: atoms += tether(m1,m1type,self.blen,self.dmin,
                           atoms[1],atoms[0],self.random)
    if m2: atoms += tether(m2,m2type,self.blen,self.dmin,
                           atoms[nsize-2],atoms[nsize-1],self.random)
    if m3: atoms += tether(m3,m3type,self.blen,self.dmin,
                           atoms[n-2],atoms[n-1],self.random)

    bonds = []
    for i in range(m1):
      if i == 0: bonds.append([1,0,n])
      else: bonds.append([1,n+i-1,n+i])
    for i in range(m2):
      if i == 0: bonds.append([1,nsize-1,n+m1])
      else: bonds.append([1,n+m1+i-1,n+m1+i])
    for i in range(m3):
      if i == 0: bonds.append([1,n-1,n+m1+m2])
      else: bonds.append([1,n+m1+m2+i-1,n+m1+m2+i])

    volume = (nsize*(nsize+1)/2 + m1+m2+m3) * pi / 6.0
    return atoms,bonds,volume

  # --------------------------------------------------------------------
  # params = m1,m2,m3,m4,m5,m6,ntype,m1type,m2type,m3type,m4type,m5type,m6type
  # ntype = type of hex center
  # m123456,m123456type = length of tethers on each corner (m = 0 = no tether)
  
  def hex(self,*params):
    m1 = params[0]
    m2 = params[1]
    m3 = params[2]
    m4 = params[3]
    m5 = params[4]
    m6 = params[5]
    ntype = params[6]
    m1type = params[7]
    m2type = params[8]
    m3type = params[9]
    m4type = params[10]
    m5type = params[11]
    m6type = params[12]
  
    atoms = []
    atoms.append([ntype,0.0,0.0,0.0])
    atoms.append([ntype,self.blen,0.0,0.0])
    atoms.append([ntype,-self.blen,0.0,0.0])
    atoms.append([ntype,self.blen/2.0,self.blen*sqrt(3.0)/2.0,0.0])
    atoms.append([ntype,-self.blen/2.0,self.blen*sqrt(3.0)/2.0,0.0])
    atoms.append([ntype,self.blen/2.0,-self.blen*sqrt(3.0)/2.0,0.0])
    atoms.append([ntype,-self.blen/2.0,-self.blen*sqrt(3.0)/2.0,0.0])
    
    n = len(atoms)
    if m1: atoms += tether(m1,m1type,self.blen,self.dmin,
                           atoms[0],atoms[1],self.random)
    if m2: atoms += tether(m2,m2type,self.blen,self.dmin,
                           atoms[0],atoms[2],self.random)
    if m3: atoms += tether(m3,m3type,self.blen,self.dmin,
                           atoms[0],atoms[3],self.random)
    if m4: atoms += tether(m4,m4type,self.blen,self.dmin,
                           atoms[0],atoms[4],self.random)
    if m5: atoms += tether(m5,m5type,self.blen,self.dmin,
                           atoms[0],atoms[5],self.random)
    if m6: atoms += tether(m6,m6type,self.blen,self.dmin,
                           atoms[0],atoms[6],self.random)

    bonds = []
    for i in range(m1):
      if i == 0: bonds.append([1,1,n])
      else: bonds.append([1,n+i-1,n+i])
    for i in range(m2):
      if i == 0: bonds.append([1,2,n+m1])
      else: bonds.append([1,n+m1+i-1,n+m1+i])
    for i in range(m3):
      if i == 0: bonds.append([1,3,n+m1+m2])
      else: bonds.append([1,n+m1+m2+i-1,n+m1+m2+i])
    for i in range(m4):
      if i == 0: bonds.append([1,4,n+m1+m2+m3])
      else: bonds.append([1,n+m1+m2+m3+i-1,n+m1+m2+m3+i])
    for i in range(m5):
      if i == 0: bonds.append([1,5,n+m1+m2+m3+m4])
      else: bonds.append([1,n+m1+m2+m3+m4+i-1,n+m1+m2+m3+m4+i])
    for i in range(m6):
      if i == 0: bonds.append([1,6,n+m1+m2+m3+m4+m5])
      else: bonds.append([1,n+m1+m2+m3+m4+m5+i-1,n+m1+m2+m3+m4+m5+i])

    volume = (7 + m1+m2+m3+m4+m5+m6) * pi / 6.0
    return atoms,bonds,volume

  # --------------------------------------------------------------------
  # params = sep,type
  # sep = separation distance of two particles
  # type = type of 2 particles

  def dimer(self,*params):
    sep = params[0]
    type = params[1]
    
    atoms = []
    x,y,z = 0.0,0.0,0.0
    atoms.append([type,x,y,z])
    x,y,z = sep,0.0,0.0
    atoms.append([type,x,y,z])

    bonds = []
    
    volume = 2 * pi / 6.0     # need to account for overlap and 2d/3d
    return atoms,bonds,volume

  # --------------------------------------------------------------------
  # params = N,sep,type
  # N = lengths of each arm of star (must be odd)
  # sep = separation distance of consecutive particles
  # type = type of all particles

  def star2d(self,*params):
    n = params[0]
    sep = params[1]
    type = params[2]
    if n % 2 == 0:
      raise StandardError, "N in patch::star2d is not odd"
    middle = n/2
    
    atoms = []
    x,y,z = 0.0,0.0,0.0
    atoms.append([type,x,y,z])
    for i in range(n):
      i -= middle
      if i == 0: continue
      x,y,z = i*sep,0.0,0.0
      atoms.append([type,x,y,z])
    for i in range(n):
      i -= middle
      if i == 0: continue
      x,y,z = 0.0,i*sep,0.0
      atoms.append([type,x,y,z])

    bonds = []
    
    volume = (2*n-1) * pi / 6.0     # need to account for overlap and 2d/3d
    return atoms,bonds,volume

  # --------------------------------------------------------------------
  # params = N,M,sep,type
  # N,M = lengths of each side of box
  # sep = separation distance of consecutive particles
  # type = type of all particles

  def box2d(self,*params):
    n = params[0]
    m = params[1]
    sep = params[2]
    type = params[3]

    height = (m-1) * sep
    width = (n-1) * sep
    
    atoms = []
    for i in range(n):
      x,y,z = i*sep,0.0,0.0
      atoms.append([type,x,y,z])
    for i in range(n):
      x,y,z = i*sep,height,0.0
      atoms.append([type,x,y,z])
    for i in range(1,m-1,1):
      x,y,z = 0.0,i*sep,0.0
      atoms.append([type,x,y,z])
    for i in range(1,m-1,1):
      x,y,z = width,i*sep,0.0
      atoms.append([type,x,y,z])

    bonds = []
    
    volume = (2*(n+m-2)) * pi / 6.0     # need to account for overlap and 2d/3d
    return atoms,bonds,volume

  # --------------------------------------------------------------------

  def random(self):
    k = self.seed/IQ
    self.seed = IA*(self.seed-k*IQ) - IR*k
    if self.seed < 0:
      self.seed += IM
    return AM*self.seed

# --------------------------------------------------------------------
# random number generator class

IM = 2147483647
AM = 1.0/IM
IA = 16807
IQ = 127773
IR = 2836

# --------------------------------------------------------------------
# push atom onto sphere surface of diam and return [type,x,y,z]

def atom_on_sphere(diam,type,x,y,z):
  if x == 0.0 and y == 0.0 and z == 0.0: scale = 0.0
  else: scale = 0.5*diam / sqrt(x*x + y*y + z*z)
  return [type,scale*x,scale*y,scale*z]

# --------------------------------------------------------------------
# make a sphere from a template with some atoms in patches w/ different types

def make_sphere(template,diam,type0,patches):
  atoms = []
  n = 0
  for atom in template:
    n += 1
    type = type0
    for pair in patches:
      patch = pair[1]
      if n in patch:
        type = pair[0]
    atoms.append(atom_on_sphere(diam,type,atom[0],atom[1],atom[2]))
  return atoms

# --------------------------------------------------------------------
# build a tether of length M of mtype, connected to atom1 with prev atom0
# blen = bond length between successive tether atoms
# dmin = length restriction between atoms i-1,i+1

def tether(m,mtype,blen,dmin,atom0,atom1,random):
  atoms = [atom0,atom1]
  for i in range(m):
    imonomer = len(atoms)
    restriction = True
    while restriction:
      rsq = 2.0
      while rsq > 1.0:
        dx = 2.0*random() - 1.0
        dy = 2.0*random() - 1.0
        dz = 2.0*random() - 1.0
        rsq = dx*dx + dy*dy + dz*dz
      r = sqrt(rsq)
      dx,dy,dz = dx/r,dy/r,dz/r
      x = atoms[imonomer-1][1] + dx*blen
      y = atoms[imonomer-1][2] + dy*blen
      z = atoms[imonomer-1][3] + dz*blen
      dx = x - atoms[imonomer-2][1]
      dy = y - atoms[imonomer-2][2]
      dz = z - atoms[imonomer-2][3]
      restriction = False
      if sqrt(dx*dx + dy*dy + dz*dz) <= dmin: restriction = True
      if not restriction: atoms.append([mtype,x,y,z])
  return atoms[2:]

# --------------------------------------------------------------------
# templates

TRI5_HOLLOW = ((0,0,0),(1,0,0),(2,0,0),(3,0,0),(4,0,0),
               (0.5,sqrt(3)/2,0),(3.5,sqrt(3)/2,0),
               (1.0,2*sqrt(3)/2,0),(3.0,2*sqrt(3)/2,0),
               (1.5,3*sqrt(3)/2,0),(2.5,3*sqrt(3)/2,0),
               (2.0,4*sqrt(3)/2,0))
               
SIMPLE_7 = ((0,0,0),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))

# C60 with added center point at end

BUCKY_60 = ((0.014416695,    1.232517,      -3.260650),
            (1.236417,       1.719517,      -2.768650),
            (2.292417,       0.8265171,     -2.492650),
            (2.122416,      -0.5504832,     -2.709650),
            (0.8954167,     -1.040483,      -3.204650),
            (-0.1565833,     -0.1504831,     -3.479650),
            (-1.018583,       1.984517,      -2.678650),
            (-0.4355836,      2.936517,      -1.826650),
            (0.9584165,      2.773517,      -1.881650),
            (1.734416,       2.937517,      -0.7156501),
            (3.065417,       0.9905167,     -1.331650),
            (2.790417,      -1.238483,      -1.682650),
            (2.234416,      -2.418483,      -1.146650),
            (1.012417,      -2.905483,      -1.638650),
            (0.3414168,     -2.215483,      -2.669650),
            (-1.052583,      -2.051483,      -2.614650),
            (-1.360583,      -0.7754831,     -3.114650),
            (-2.397583,      -0.020483017,   -2.529650),
            (-2.227583,       1.356517,      -2.312650),
            (-1.059583,       3.265517,      -0.6046500),
            (-2.262583,       2.640517,      -0.2406499),
            (-2.848583,       1.684517,      -1.095650),
            (-3.402583,       0.5095167,     -0.5616500),
            (-3.124583,      -0.5444832,     -1.447650),
            (-2.815583,      -1.825483,      -0.9456501),
            (-1.781583,      -2.577483,      -1.527650),
            (-1.113584,      -3.265483,      -0.5006499),
            (0.2864165,     -3.429483,      -0.5566499),
            (3.373417,      -0.2854834,     -0.8306501),
            (1.018416,      -1.984483,       2.678350),
            (2.227417,      -1.356483,       2.312350),
            (2.398417,       0.020516872,    2.530350),
            (1.360417,       0.7755170,      3.114350),
            (0.1564164,      0.1505170,      3.479350),
            (-1.236583,      -1.719483,       2.768350),
            (-0.9585834,     -2.773483,       1.882350),
            (0.4354167,     -2.936483,       1.826350),
            (1.059417,      -3.265483,       0.6053500),
            (2.263417,      -2.640483,       0.2403500),
            (2.848417,      -1.684483,       1.096350),
            (3.124417,       0.5445170,      1.447350),
            (2.815417,       1.825517,       0.9453500),
            (1.781417,       2.577517,       1.527350),
            (1.052417,       2.051517,       2.614350),
            (-0.3415833,      2.215517,       2.669350),
            (-0.8955836,      1.040517,       3.204350),
            (-2.122583,       0.5505171,      2.710350),
            (-2.292583,      -0.8264832,      2.492350),
            (-1.734583,      -2.937483,       0.7163500),
            (-2.786583,      -2.048483,       0.4413500),
            (-3.065583,      -0.9904833,      1.331350),
            (-3.373584,       0.2855167,      0.8313500),
            (-2.790584,       1.238517,       1.683350),
            (-2.233583,       2.418517,       1.146350),
            (-1.011583,       2.905517,       1.638350),
            (-0.2855835,      3.429517,       0.5563500),
            (1.113417,       3.265517,       0.5003500),
            (3.402417,      -0.5094833,      0.5613500),
            (2.786417,       2.047517,      -0.4416499),
            (-0.014583588,   -1.232483,       3.261350),
            (0.0,            0.0,             0.0))

# --------------------------------------------------------------------
# C80

BUCKY_80 = ((-1.243762,     -0.7016125,       3.734700),
            (-1.248762,      0.6973875,       3.749700),
            (-0.054762483,   1.429388,       3.749700),
            (1.189238,      0.7893875,       3.735700),
            (1.194237,     -0.6086125,       3.783700),
            (0.00023752451,  -1.340613,       3.782700),
            (-2.204762,       1.428387,       3.034700),
            (-3.191762,      0.7883875,       2.276700),
            (-3.223763,     -0.6096125,       2.311700),
            (-2.267762,      -1.340613,       3.027700),
            (0.2282375,      -2.540613,       3.098700),
            (-0.7787625,      -3.147613,       2.340700),
            (-2.039762,      -2.541613,       2.343700),
            (2.160238,      -1.356613,       3.098700),
            (1.563237,      -2.550612,       2.675700),
            (-0.2727625,       2.612388,       3.034700),
            (-1.601763,       2.612388,       2.591700),
            (-3.586762,      -1.357612,       1.184700),
            (-2.854763,      -2.551612,       1.204700),
            (2.223238,       1.412387,       3.028700),
            (2.006238,       2.595387,       2.312700),
            (0.7452375,       3.202387,       2.276700),
            (0.3742375,       3.844388,       1.090700),
            (-0.9547625,       3.844388,      0.6487002),
            (-1.962763,       3.201387,       1.374700),
            (-2.992763,       2.595387,      0.6477003),
            (-3.595763,       1.411387,       1.090700),
            (-3.958763,      0.6633875,     -0.036299706),
            (-3.930763,     -0.7356125,     -0.020299911),
            (-3.568762,      -1.385612,      -1.205300),
            (-2.837763,      -2.579613,      -1.186300),
            (-2.439763,      -3.168612,      0.019700050),
            (-1.205762,      -3.828612,      0.035700321),
            (-0.3907624,      -3.817612,       1.174700),
            (0.9442375,      -3.827612,      0.7516999),
            (1.941237,      -3.168612,       1.478700),
            (3.190238,      0.6643875,       2.343700),
            (3.158237,     -0.7346125,       2.340700),
            (3.579237,      -1.384613,       1.175700),
            (2.982238,      -2.578612,      0.7527003),
            (1.244238,      0.7023875,      -3.735300),
            (1.249238,     -0.6976125,      -3.750300),
            (0.055237532,  -1.428612,      -3.750300),
            (-1.189762,     -0.7896125,      -3.735300),
            (-1.193763,      0.6093875,      -3.783300),
            (0.00023752451,   1.340387,      -3.783300),
            (2.205238,      -1.428612,      -3.034300),
            (3.192237,     -0.7886125,      -2.276300),
            (3.224237,      0.6093875,      -2.312300),
            (2.268238,       1.341388,      -3.027300),
            (-0.2277625,       2.541388,      -3.098300),
            (0.7792375,       3.147388,      -2.340300),
            (2.040237,       2.541388,      -2.343300),
            (-2.159763,       1.356387,      -3.098300),
            (-1.562762,       2.550387,      -2.675300),
            (0.2732375,      -2.612613,      -3.034300),
            (1.601238,      -2.612613,      -2.592300),
            (3.586237,       1.357388,      -1.185300),
            (2.854238,       2.551387,      -1.204300),
            (-2.223763,      -1.411613,      -3.028300),
            (-2.005763,      -2.595613,      -2.312300),
            (-0.7457625,      -3.201612,      -2.277300),
            (-0.3747625,      -3.844613,      -1.090300),
            (0.9542375,      -3.844613,     -0.6482999),
            (1.962238,      -3.201612,      -1.375300),
            (2.992238,      -2.595613,     -0.6482999),
            (3.596237,      -1.411613,      -1.090300),
            (3.958237,     -0.6636125,      0.036700249),
            (3.931237,      0.7353874,      0.019700050),
            (3.569237,       1.385387,       1.205700),
            (2.837237,       2.579388,       1.185700),
            (2.439238,       3.168387,     -0.019299984),
            (1.205238,       3.828387,     -0.036299706),
            (0.3902375,       3.818388,      -1.175300),
            (-0.9447625,       3.828387,     -0.7522998),
            (-1.941762,       3.168387,      -1.478300),
            (-3.189763,     -0.6646125,      -2.344300),
            (-3.157763,      0.7343875,      -2.340300),
            (-3.579762,       1.384387,      -1.175300),
            (-2.982763,       2.578387,     -0.7522998))
