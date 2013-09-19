""" reaxFF ffield file format parser
@status: Read the reactive force field (field) and catalog the parameters 
         into functions as JPCB 2008.
         05-30: support Inner-wall and lg
@todo:
To output the ffield
@see: output_ff.py
@note: for ATOM n index is exactly n in table

1. Combination rules:
    comb1: (a + b) / 2
    comb2: math.sqrt(a*b)
    comb3: 2*math.sqrt(a*b)
    comb4: math.pow(a*b, 1.5)
"""
import os
import math
import argparse
DEBUG = 0

class Ffield():
    """ reaxFF force field (ffield) file
    """
    def __init__(self, filename="ffield", lg=0):
        self.lg = lg
        """@ivar: lg type ffield or not"""
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
        self.elements= []
        """@ivar: elements in the ffield"""
        self.read(filename)
        """@xxxxx: read the input file """

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
            if self.lg == 0:
                param = self.readParams(4, f)
                self.atom.append(param)
            elif self.lg == 1:
                param = self.readParams(5, f)
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
        # Get elements map
        self.getMap()

    def ati(self, line):
        return int(line.strip().split()[0])

    def atf(self, line):
        return float(line.strip().split()[0])

    def getMap(self,):
        if len(self.atom) > 0:
            for i in self.atom:
                self.elements.append(i[0])

    def readParams(self, n, f):
        params = []
        for i in range(n):
            tokens = f.readline().strip().split()
            for j in tokens:
                params.append(j)
        return params

    def toEquation(self,):
        """output the force field parameter to more readable form.
        according to reaxFF 2008
        """

        # Get the bond order related parameters
        # get the r_sigma (r_s), r_pi and r_pipi using comb1
        # @note: only work for both of the atomic r > 0
        eq2 = []
        for i in range(len(self.atom)):
            for j in range(i, len(self.atom) - 1):
                r_s = 0.5 * (float(self.atom[i][1]) + float(self.atom[j][1]))
                if float(self.atom[i][7]) > 0 and float(self.atom[j][7]) > 0:
                    r_pi = 0.5 * (float(self.atom[i][7]) + float(self.atom[j][7]))
                else:
                    r_pi = -1.0
                if float(self.atom[i][17]) > 0 and float(self.atom[j][17]) > 0:
                    r_pipi = 0.5 * (float(self.atom[i][17]) + float(self.atom[j][17]))
                else:
                    r_pipi = -1.0
                eq2.append([i, j, r_s, 0, 0, r_pi, 0, 0, r_pipi, 0, 0])

        # get the p_bo(n) ( n = 1, 2, ..6) from the bond section
        # @note: only work for both of the atomic r > 0
        n, n1 = 0, 0
        for i in range(len(self.atom)):
            for j in range(i, len(self.atom) - 1):
                p_bo1, p_bo2, p_bo3 = 0, 0, 0
                p_bo4, p_bo5, p_bo6 = 0, 0, 0
                a1 = i + 1
                a2 = j + 1
                for k in self.bond:
                    if int(k[0]) == a1 and int(k[1]) == a2:
                        p_bo1, p_bo2, p_bo3 = float(k[14]), float(k[15]), float(k[11])
                        p_bo4, p_bo5, p_bo6 = float(k[12]), float(k[6]), float(k[8])
                        n1 += 1
                    elif int(k[0]) == a2 and int(k[1]) == a1:
                        p_bo1, pbo2, pbo3 = float(k[14]), float(k[15]), float(k[11])
                        p_bo4, pbo5, pbo6 = float(k[12]), float(k[6]), float(k[8])
                        n1 += 1
                eq2[n][3], eq2[n][4] = p_bo1, p_bo2
                eq2[n][6], eq2[n][7] = p_bo3, p_bo4
                eq2[n][9], eq2[n][10] = p_bo5, p_bo6
                n += 1

        # update the r_sigma, r_pi and r_pipi according to the off-diagonal parameters
        n, n1 = 0, 0
        for i in range(len(self.atom)):
            for j in range(i, len(self.atom) - 1):
                a1 = i + 1
                a2 = j + 1
                for k in self.off:
                    if int(k[0]) == a1 and int(k[1]) == a2:
                        eq2[n][2] = float(k[5])
                        eq2[n][5] = float(k[6])
                        eq2[n][8] = float(k[7])
                        n1 += 1
                    elif int(k[0]) == a2 and int(k[1]) == a1:
                        eq2[n][2] = float(k[5])
                        eq2[n][5] = float(k[6])
                        eq2[n][8] = float(k[7])
                        n1 += 1
                n += 1

        # equation 4
        eq4 = []
        for i in range(len(self.atom)):
            for j in range(i, len(self.atom) - 1):
                p_o_131 = math.sqrt(float(self.atom[i][20])*float(self.atom[j][20]))
                p_o_132 = math.sqrt(float(self.atom[i][21])*float(self.atom[j][21]))
                p_o_133 = math.sqrt(float(self.atom[i][22])*float(self.atom[j][22]))
                eq4.append([i, j, p_o_132, p_o_131, p_o_133])
        # equation 6
        eq6 = []
        n1 = 0
        for i in range(len(self.atom)):
            for j in range(i, len(self.atom) - 1):
                de_s, de_p, de_pp = 0, 0, 0
                p_be1, p_be2 = 0, 0
                a1 = i + 1
                a2 = j + 1
                for k in self.bond:
                    if int(k[0]) == a1 and int(k[1]) == a2:
                        de_s, de_p, de_pp = float(k[2]), float(k[3]), float(k[4])
                        p_be1, p_be2 = float(k[5]), float(k[10])
                        n1 += 1
                    elif int(k[0]) == a2 and int(k[1]) == a1:
                        de_s, de_p, de_pp = float(k[2]), float(k[3]), float(k[4])
                        p_be1, p_be2 = float(k[5]), float(k[10])
                        n1 += 1
                eq6.append([i, j, de_s, p_be1, p_be2, de_p, de_pp])

        # equation 23
        eq23 = []
        for i in range(len(self.atom)):
            for j in range(i, len(self.atom) - 1):
                epsilon = math.sqrt(float(self.atom[i][5])*float(self.atom[j][5]))
                alpha = math.sqrt(float(self.atom[i][9])*float(self.atom[j][9]))
                #r_vdw = 2 * math.sqrt(float(self.atom[i][4])*float(self.atom[j][4]))
                r_vdw = math.sqrt(float(self.atom[i][4])*float(self.atom[j][4]))
                gamma_w = 0.5 * (float(self.atom[i][10]) + float(self.atom[j][10]))
                eq23.append([i, j, epsilon, alpha, r_vdw, gamma_w])

        n, n1 = 0, 0
        for i in range(len(self.atom)):
            for j in range(i, len(self.atom) - 1):
                a1 = i + 1
                a2 = j + 1
                for k in self.off:
                    if int(k[0]) == a1 and int(k[1]) == a2:
                        eq23[n][2] = float(k[2])
                        eq23[n][3] = float(k[4])
                        eq23[n][4] = float(k[3])
                        n1 += 1
                    elif int(k[0]) == a2 and int(k[1]) == a1:
                        eq23[n][2] = float(k[2])
                        eq23[n][3] = float(k[4])
                        eq23[n][4] = float(k[3])
                        n1 += 1
                n += 1

        # equation 26
        eq26 = []
        for i in range(len(self.atom)):
            for j in range(i, len(self.atom) - 1):
                d = math.sqrt(float(self.atom[i][31])*float(self.atom[j][31]))
                #print self.atom[i][0], self.atom[j][0], self.atom[i][31], self.atom[j][31], d
                alpha = math.sqrt(float(self.atom[i][32])*float(self.atom[j][32]))
                r_inner = math.sqrt(float(self.atom[i][30])*float(self.atom[j][30]))
                eq26.append([i, j, d, alpha, r_inner])

        # output
        o = open("output.ff", "w")
        o.write("# equation 2: bond order\n")
        o.write("%4s%4s%8s%8s%8s"%("a1", "a2", "ro_s", "pbo1", "pbo2"))
        o.write("%8s%8s%8s%8s%8s%8s\n"%("ro_pi", "pbo3", "pbo4", "ro_pipi", "pbo5", "pbo6"))
        for i in range(len(eq2)):
            counter = 0 
            for j in range(len(eq2[i])):
                if counter < 2:
                    o.write("%4s"%self.elements[eq2[i][j]])
                else:
                    o.write("%8.4f"%eq2[i][j])
                counter += 1
            o.write("\n")
        
        o.write("# equation 4: bond order corrections\n")
        o.write("%4s%4s%8s%8s%8s\n"%("a1", "a2", "pboc3", "pboc4", "pboc5"))
        for i in range(len(eq4)):
            counter = 0 
            for j in range(len(eq4[i])):
                if counter < 2:
                    o.write("%4s"%self.elements[eq4[i][j]])
                else:
                    o.write("%8.4f"%eq4[i][j])
                counter += 1
            o.write("\n")

        o.write("# equation 6: bond energy\n")
        o.write("%4s%4s%10s%10s"%("a1", "a2", "de_s", "p_be1"))
        o.write("%10s%10s%10s\n"%("p_be2", "de_p", "de_pp"))
        for i in range(len(eq6)):
            counter = 0 
            for j in range(len(eq6[i])):
                if counter < 2:
                    o.write("%4s"%self.elements[eq6[i][j]])
                else:
                    o.write("%10.4f"%eq6[i][j])
                counter += 1
            o.write("\n")

        o.write("# equation 23: van der Waals interaction\n")
        o.write("%4s%4s%10s%10s"%("a1", "a2", "epsil", "alpha"))
        o.write("%10s%10s\n"%("r_vdw", "gamma_w"))
        for i in range(len(eq23)):
            counter = 0 
            for j in range(len(eq23[i])):
                if counter < 2:
                    o.write("%4s"%self.elements[eq23[i][j]])
                else:
                    o.write("%10.4f"%eq23[i][j])
                counter += 1
            o.write("\n")

        """
        o.write("# equation 24: coulomb interaction\n")
        o.write("%4s%4s%10s%10s"%("a1", "a2", "chiEEM", "etaEEM"))
        o.write("%10s\n"%("gamma_w"))
        """

        o.write("# equation 26: inner wall\n")
        o.write("%4s%4s%10s%10s"%("a1", "a2", "d", "alpha"))
        o.write("%10s\n"%("r_inner"))
        for i in range(len(eq26)):
            counter = 0 
            for j in range(len(eq26[i])):
                if counter < 2:
                    o.write("%4s"%self.elements[eq26[i][j]])
                else:
                    o.write("%10.4f"%eq26[i][j])
                counter += 1
            o.write("\n")

        o.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="ffield", nargs="?", help="force field file name")
    parser.add_argument("-D", action="store_true", help="Debug the code")
    parser.add_argument("-type", nargs=1, type=int, help="Force field type: 0 for vdw; 1 for lg_inner wall")
    args = parser.parse_args()
    #print b.getBondDist(3,2)
    
    fname = args.fname

    assert os.path.exists(fname)

    if args.D:
        DEBUG = 1

    if args.type:
        ntype = args.type[0]
    else:
        ntype = 0
        print "Warning: Using default force field type"

    ff = Ffield(fname, ntype)
    ff.toEquation()

if __name__ == "__main__":
    main()
