"""basic data structures of molecule and atoms
@version: 1.0
@author: hawkweedcheng
@contact: chengtao@sjtu.edu.cn
@todo: 
        - 2016-10-05 convert xyz to frac;
"""
import numpy
import math
import re
from cons import ELEMENT2MASS 
from utilities import get_dist, get_angle

class Atom():
    """Basic class for atom 
    @todo: extend to general class
    """
    def __init__(self, ):
        self.number = 0
        """@ivar: atom number from input file
        @type: int
        """
        self.name = ''
        """@ivar: atom name from input file
        @type: char
        """
        self.an = 0
        """@ivar: atomic number
        @type: int
        """
        self.type1 = 0
        """@ivar: apparent atom type
        @type: int
        """
        self.type2 = 0
        """@ivar: default atom type (DFF)
        @type: int
        """
        self.resn = 0
        """@ivar: residue number belong to
        @type: char
        """
        self.resname = ''
        """@ivar: resname belong to
        @type: char
        """
        self.x = [0, 0, 0]
        """@ivar: atom coordination
        @type: list
        """
        self.xFrac = [0, 0, 0]
        """@ivar: atom fraction coordination
        @type: list
        """
        self.xr = [0, 0, 0]
        """@ivar: restraint in x, y and z direction
        @type: list
        """
        self.element = ''
        """@ivar: element
        @type: char 
        """
        self.charge = 0.0
        """@ivar: charge
        @type: float
        """
        self.cg = 0
        """@ivar: charge group belong to
        @type: int
        """

class Bond():
    """Basic class for bond
    """
    def __init__(self,):
        self.b1 = 0
        """@ivar: atom number of atom1 in bond"""
        self.b2 = 0
        """@ivar: atom number of atom2 in bond"""
        self.type = 0
        """@ivar: bond type number"""
    def getDistance(self,):
        """Return bond distance
        @todo: doing.....
        """
        pass

class System():
    """Basic class for system
    @todo:
    """
    def __init__(self,):
        self.name = ''
        """@ivar: system name
        @type: char
        """
        self.step = 0

        self.atoms = []
        self.natoms = []
        self.n_atoms = 0
        self.atomtypes = []
        self.n_atomtypes = 0
        self.bonds = []
        self.n_bonds = 0
        self.n_bondtypes = 0
        self.angles = []
        self.n_angles = 0
        self.n_angletypes = 0
        self.dihedrals = []
        self.n_dihedrals = 0
        self.n_dihedraltypes= 0
        self.impropers= []
        self.n_impropers= 0
        self.n_impropertypes= 0
        """@ivar: atom types in system
        @type: list
        """
        self.mols = []
        """@ivar: molecules in system
        @type: list
        @note: parseToMol()"""
        self.mass = 0.0
        """@ivar: total atomic mass of system
        @type: float
        """
        self.mlist = {}
        """@ivar: molecule list (format {resname: [ 1, 2, ....]})
        @type: dict 
        @note: """
        self.pbc = []
        """@ivar: pbc 
        @type: list 
        @note: """
        self.scaleFactor = 1.0
        """@ivar: scale factor
        @type: float
        @note: """
        self.map = []
        """@ivar: the element-atom type map
        @type: list 
        @note: """
        #specify for geo file
        self.geotag = ''
        """@ivar: maps
        @type: list 
        @note: """
        #specify for g03 gjf
        self.options = []
        """@ivar: options in g03, such as memory
        @type: list 
        @note: """
        self.methods = []
        """@ivar: computational methods, such as opt
        @type: list 
        @note: """
        self.connect = []
        """@ivar: connectivity in gjf file
        @type: list 
        @note: should make a standary here
        """
        self.redundant = []
        """@ivar: redundant info, such as fix bond distances
        @type: list 
        @note: """
        self.spin = 1
        """@ivar: spin
        @type: int
        @note: """
        self.charge = 0
        """@ivar: charge
        @type: int
        @note: """
        self.ffparams = []
        """@ivar: force field parameters
        @type: list
        """

    def parseToMol(self,):
        """Catalog atoms in system to molecules according to resname. 
        The bonds also parsed into each molecules"""

        #catalog atoms according to resname
        for i in self.atoms:
            if i.resn in self.mlist.keys():
                m.atoms.append(i)
                self.mlist[i.resn].append(i.number)
            else:
                if len(self.mlist.keys()) > 0:
                    self.mols.append(m)
                m = Molecule()
                m.name = i.resn
                m.atoms.append(i)
                self.mlist[i.resn] = []
                self.mlist[i.resn].append(i.number)
        self.mols.append(m)
        self.mlist[i.resn].append(i.number)

        #catalog bonds according to molecule list (self.mlist)
        for i in self.mols:
            for j in self.bonds:
                if j.b1 in self.mlist[i.name] and j.b2 in self.mlist[i.name]:
                    i.bonds.append(j)

    def assignEleTypes(self,):
        """ assign element types 
        """
        pattern = re.compile(r'(\D+)(\d*)')
        counter = 0
        for i in self.atoms:
            b = i.name.strip()
            if "_" in b:
                b = b.split("_")[0]
            match = pattern.match(b)
            if match:
                i.element= match.group(1)
                if len(i.element) > 1:
                    i.element = i.element[0].upper() + i.element[1].lower()
    
    def assignAtomTypes2(self,):
        """ assign atomtypes according to the element types
        specially designed for lammps data file
        """
        self.assignEleTypes()
        for i in self.atoms:
            i.name = i.element

    def assignAtomTypes(self,):
        """ assign atomtypes according to the element types
        specially designed for lammps data file
        """
        a = []
        c = {}
        counter = 0
        for i in self.atoms:
            b = i.name.strip()
            if b not in a:
                counter += 1
                a.append(b)
                c[b] = counter
                self.map.append([counter, b])
        for i in self.atoms:
            i.type1 = c[i.name.strip()]

    def assignIdNumbers(self,):
        """ assign id numbers according to the sequence of atoms
            starting from 1
        """
        n = 1
        for i in self.atoms:
            i.number = n
            n += 1

    def toFrac(self,):
        a2r = numpy.pi / 180.0
        alpha = a2r * self.pbc[3]
        beta = a2r * self.pbc[4]
        gamma = a2r * self.pbc[5]
        a, b, c = self.pbc[0], self.pbc[1], self.pbc[2]
        v = math.sqrt(1 -math.cos(alpha)*math.cos(alpha)\
                 - math.cos(beta)*math.cos(beta)\
                 - math.cos(gamma)*math.cos(gamma)\
                 + 2*math.cos(alpha)*math.cos(beta)*math.cos(gamma))
        tmat = numpy.matrix( [
            [1.0/a, -math.cos(gamma)/(a*math.sin(gamma)), (math.cos(alpha)*math.cos(gamma)-math.cos(beta))/(a*v*math.sin(gamma))],
            [ 0.0, 1.0/(b*math.sin(gamma)), (math.cos(beta)*math.cos(gamma)-math.cos(alpha))/ (b*v*math.sin(gamma))  ],
            [ 0.0, 0.0, math.sin(gamma)/(c*v)]])
        for i in self.atoms:
            cart_coord = numpy.array(i.x)
            r = cart_coord*tmat.T
            i.xFrac = numpy.squeeze(numpy.asarray(r))
        
    def translate(self, delta=0.0, axis="z"):
        """ translate the coordination along x, y or z
        """
        if axis == "x":
            n = 0
        elif axis == "y":
            n = 1
        elif axis == "z":
            n = 2
        for i in self.atoms:
            i.x[n] = i.x[n] + delta
    
    def scale(self, scale=1.0):
        """ scale the  coordinations
        """
        for i in self.atoms:
            i.x[n] = i.x[n] * scale

    def getMin(self, axis="z"):
        """ get the min value in x, y or z
        """
        min = 9999.0
        if axis == "x":
            n = 0
        elif axis == "y":
            n = 1
        elif axis == "z":
            n = 2
        for i in self.atoms:
            val = i.x[n]
            if  min > i.x[n]:
                min = i.x[n]
        return min

    def getVol(self):
        """ return the volume of the system
        @see: http://en.wikipedia.org/wiki/Parallelepiped
        """
        vol = 0
        if len(self.pbc) > 0:
            a = self.pbc[0]
            b = self.pbc[1]
            c = self.pbc[2]
            alpha = self.pbc[3]/180.0*math.pi
            beta = self.pbc[4]/180.0*math.pi
            gamma = self.pbc[5]/180.0*math.pi
            cos = 2*math.cos(alpha)*math.cos(beta)*math.cos(gamma)
            cos2 = math.cos(alpha)*math.cos(alpha) + math.cos(beta)*math.cos(beta) \
                    + math.cos(gamma)*math.cos(gamma)
            sq = math.sqrt(1+2*cos-cos2)
            vol = a*b*c*sq
        else:
            print("Warning: no box defined")
        return vol

    def getMass(self):
        """Get the mass
        """
        self.assignEleTypes()
        mass = 0.0
        for i in self.atoms:
            mass += ELEMENT2MASS[i.element] 
        self.mass = mass
        return mass

    def getBondDist(self, a1, a2):
        """Get the distance
        """
        x1 = self.atoms[a1-1].x
        x2 = self.atoms[a2-1].x
        dist = get_dist(x1, x2)
        return dist

    def getAngle(self, a1, a2, a3):
        """Get the angle
        """
        x1 = self.atoms[a1-1].x
        x2 = self.atoms[a2-1].x
        x3 = self.atoms[a3-1].x
        ang = get_angle(x1, x2, x3)
        return ang
        
    def sortNdx(self, ndx):
        atoms = []
        for i in ndx:
            atoms.append(self.atoms[i])
        self.atoms = atoms

    def sortXYZ(self, axis="z"):
        atoms = []
        if axis == "x":
            atoms = sorted(self.atoms, key=lambda atom: atom.x[0])
        elif axis == "y":
            atoms = sorted(self.atoms, key=lambda atom: atom.x[1])
        elif axis == "z":
            atoms = sorted(self.atoms, key=lambda atom: atom.x[2])
        else:
            print("Error: non-sense input")
            atoms = self.atoms

        self.atoms = atoms

class Molecule():
    """Basic class for molecular includes (name, atoms(list of atoms), 
    conn (connectivity))
    @todo:
    """
    def __init__(self, ):
        self.name = ''
        self.n_atoms = 0
        self.n_bonds = 0
        self.n_angles = 0
        self.n_dihedrals = 0
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.a0 = 0
        self.a1 = 0

        self.nb = {}
        """@ivar: neighbour list
        @type: dict
        """

    def getNB(self):
        """Get neighbour list from bond list according to atom number
        @todo:
        """
        if self.bonds:
            for i in self.atoms:
                n = i.number
                if n not in self.nb.keys():
                    self.nb[n] = []
                for j in self.bonds:
                    if j.b1 == n:
                        self.nb[n].append(j.b2)
                    elif j.b2 == n:
                        self.nb[n].append(j.b1)
                    else:
                        pass
        else:
            print("Note : There is no bond info in inputfile!")

    def getMaxCoord(self, n=0):
        """return the maximun coordinate value in x, y or z ( n = 0, 1, 2)
        direction. Default is x direction"""
        max = -10000
        for i in self.mol:
            if max < i.x[n]:
                max = i.x[n]
        print(max)

    def getMinCoord(self, n=0):
        """return the min coordinate value in x, y or z ( n = 0, 1, 2)
        direction. Default is x direction"""
        min = 10000
        for i in self.mol:
            if min > i.x[n]:
                min = i.x[n]
        print(min)
    
def test():
    """test some functions
    """ 
    pass

if __name__ == "__main__":
    test()

