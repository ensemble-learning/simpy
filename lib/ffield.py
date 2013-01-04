""" reaxFF ffield file format parser
"""
DEBUG = 0

class Ffield():
    """ reaxFF force field (ffield) file
    """
    def __init__(self, filename="ffield"):
        self.gl = []
        """@ivar: global parameters"""
        self.atom = []
        """@ivar: atom parameters"""
        self.bond = []
        """@ivar: bond parameters"""
        self.off = []
        """@ivar: off parameters"""
        self.angle = []
        """@ivar: angle parameters"""
        self.torsion = []
        """@ivar: torsion parameters"""
        self.hbond = []
        """@ivar: hbond parameters"""
        self.read(filename)

    def read(self, filename):
        """ read the ffield file
        """
        f = open(filename, 'r')

        # read the comment
        f.readline()

        # GLOBAL parameters
        if DEBUG:
            print "-"*20, "reading GLOBAL parameters", "-"*20
        n = self.ati(f.readline())
        self.gl = []
        for i in range(n):
            self.gl.append(self.atf(f.readline()))

        # ATOM parameters
        if DEBUG:
            print "-"*20, "reading ATOM parameters", "-"*20
        n = self.ati(f.readline())    
        f.readline()
        f.readline()
        f.readline()
        self.atom = []
        for i in range(n):
            param = self.readParams(4, f)
            self.atom.append(param)

        # BOND parameters
        if DEBUG:
            print "-"*20, "reading BOND parameters", "-"*20
        n = self.ati(f.readline())    
        f.readline()
        self.bond = []
        for i in range(n):
            param = self.readParams(2, f)
            self.bond.append(param)

        # OFF parameters
        if DEBUG:
            print "-"*20, "reading OFF parameters", "-"*20
        n = self.ati(f.readline())    
        self.off = []
        for i in range(n):
            param = self.readParams(1, f)
            self.off.append(param)

        # ANGLE parameters
        if DEBUG:
            print "-"*20, "reading ANGLE parameters", "-"*20
        n = self.ati(f.readline())    
        self.angle = []
        for i in range(n):
            param = self.readParams(1, f)
            self.angle.append(param)

        # TORSION parameters
        if DEBUG:
            print "-"*20, "reading TORSION parameters", "-"*20
        n = self.ati(f.readline())    
        self.torsion = []
        for i in range(n):
            param = self.readParams(1, f)
            self.torsion.append(param)

        # H-BOND parameters
        if DEBUG:
            print "-"*20, "reading H-BOND parameters", "-"*20
        n = self.ati(f.readline())    
        self.hbond = []
        for i in range(n):
            param = self.readParams(1, f)
            self.hbond.append(param)
        f.close()

    def ati(self, line):
        return int(line.strip().split()[0])

    def atf(self, line):
        return float(line.strip().split()[0])

    def readParams(sefl, n, f):
        params = []
        for i in range(n):
            tokens = f.readline().strip().split()
            for j in tokens:
                params.append(j)
        return params

if __name__ == "__main__":
    ff = Ffield("ffield2")
    print ff.atom
    
