"""utilities
"""
import math

def t2element():
    pass


def rules(line, keyword):
    if len(line.strip()) == 0:
        return 1
    else:
        return 0

def parseBlock(filename, keyword):
    """ a general parse into blocks
    """
    flag = 1
    blocks = []
    block = []
    f = open(filename, "r")
    
    for i in f:
        flag = 1
        if rules(i, keyword):
            if len(block) > 0:
                block.append(i)
                blocks.append(block)
                block = []
        else:
            block.append(i)

    for i in range(len(blocks)):
        o = open("out%03d"%i, 'w')
        for j in blocks[i]:
            o.write(j)
        o.close()

def toDegree(a):
    return a/math.pi*180.0

def v2lattice(x, y, z):
    """convert vectors to a, b, c, alpha, beta and gamma
    @see: Triclinic (non-orthogonal) simulation boxes (http://lammps.sandia.gov/doc/Section_howto.html)
    """
    a = x[0]
    b = math.sqrt(y[1]*y[1]+y[0]*y[0])
    c = math.sqrt(z[2]*z[2] + z[0]*z[0] + z[1]*z[1])
    alpha = math.acos((y[0]*z[0] + y[1]*z[1])/(b*c))
    alpha = toDegree(alpha)
    beta = math.acos(z[0]/c)
    beta = toDegree(beta)
    gamma = math.acos(y[0]/b)
    gamma = toDegree(gamma)
    return [a, b, c, alpha, beta, gamma]

def lattice2v(pbc):
    """convert a, b, c, alpha, beta and gamma to vectors
    @see: v2lattice()
    """
    a = pbc[0]
    b = pbc[1]
    c = pbc[2]
    alpha = pbc[3] / 180.0*math.pi
    beta = pbc[4] / 180.0*math.pi
    gamma = pbc[5] / 180.0*math.pi
    xx = a
    xy = b * math.cos(gamma)
    xz = c * math.cos(beta)
    yy = math.sqrt(b * b - xy * xy)
    yz = (b * c * math.cos(alpha) - xy * xz) / yy
    zz = math.sqrt(c * c - xz * xz - yz * yz)
    
    return xx, xy, xz, yy, yz, zz

def get_dist(x1, x2):
    a = x1[0] - x2[0]
    b = x1[1] - x2[1]
    c = x1[2] - x2[2]
    return math.sqrt(a*a + b*b + c*c)

if __name__ == "__main__":
    a =  [3.4199999522999995 ,  0.0000000000000000  , 0.0000000000000000]
    b =  [-1.7099999761999998 ,  2.9608964930999999 ,  0.0000000000000000]
    c =  [0.0000000000000000  , 0.0000000000000000  , 5.6100000191000001]
    a =  [25.50 ,  0.0000000000000000  , 0.0000000000000000]
    b =  [-12.75,  21.91,  0.0000000000000000]
    c =  [0.0000000000000000  , 0.0000000000000000  , 24.88]
    pbc = v2lattice(a, b, c)
    #print pbc

