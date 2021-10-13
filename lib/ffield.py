""" reaxFF ffield file format parser
@status: Read the reactive force field (field) and catalog the parameters 
         into functions as JPCB 2008.
         05-30: support Inner-wall and lg
@todo:
@0108_2014: to complete the checkRedudant()
@Sun Mar 23 23:06:46 PDT 2014
    - check out & merge

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
from output_ff import toFfield
DEBUG = 1

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
            print("-"*20, "reading GLOBAL parameters", "-"*20)
        n = self.ati(f.readline())
        self.gl = []
        for i in range(n):
            self.gl.append(self.atf(f.readline()))

        # ATOM parameters
        if DEBUG:
            print("-"*20, "reading ATOM parameters", "-"*20)
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
            print("-"*20, "reading BOND parameters", "-"*20)
        n = self.ati(f.readline())    
        f.readline()
        self.bond = []
        for i in range(n):
            param = self.readParams(2, f)
            self.bond.append(param)

        # OFF parameters
        if DEBUG:
            print("-"*20, "reading OFF parameters", "-"*20)
        n = self.ati(f.readline())    
        self.off = []
        for i in range(n):
            param = self.readParams(1, f)
            self.off.append(param)

        # ANGLE parameters
        if DEBUG:
            print("-"*20, "reading ANGLE parameters", "-"*20)
        n = self.ati(f.readline())    
        self.angle = []
        for i in range(n):
            param = self.readParams(1, f)
            self.angle.append(param)

        # TORSION parameters
        if DEBUG:
            print("-"*20, "reading TORSION parameters", "-"*20)
        n = self.ati(f.readline())    
        self.torsion = []
        for i in range(n):
            param = self.readParams(1, f)
            self.torsion.append(param)

        # H-BOND parameters
        if DEBUG:
            print("-"*20, "reading H-BOND parameters", "-"*20)
        n = self.ati(f.readline())    
        self.hbond = []
        for i in range(n):
            param = self.readParams(1, f)
            self.hbond.append(param)
        f.close()
        # Get elements map
        self.getMap()
    
    def clearup(self,):
        self.atom = []
        self.bond = []
        self.off = []
        self.angle = []
        self.torsion = []
        self.hbond = []
        self.elements= []

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

    def completeOff(self,):
        """ Complete the off-dia section using combination rule
        """
        off_ext = []
        for i in range(len(self.atom)):
            for j in range(i+1, len(self.atom)):
                a1 = i + 1
                a2 = j + 1
                flag = 0
                for k in self.off:
                    off1 = int(k[0])
                    off2 = int(k[1])
                    if (a1 == off1 and a2 == off2) or (a1 == off2 and a2 == off1):
                        flag += 1
                if flag == 0:
                    atom1 = str(a1)
                    atom2 = str(a2)
                    Dij = math.sqrt(float(self.atom[i][5])*float(self.atom[j][5]))
                    RvdW = math.sqrt(float(self.atom[i][4])*float(self.atom[j][4]))
                    alpha = math.sqrt(float(self.atom[i][9])*float(self.atom[j][9]))
                    r_s = (float(self.atom[i][1])*float(self.atom[j][1]))/2.0
                    r_pi = (float(self.atom[i][7])*float(self.atom[j][7]))/2.0
                    r_pipi = (float(self.atom[i][17])*float(self.atom[j][17]))/2.0
                    off_ext.append([atom1, atom2, Dij, RvdW, alpha, r_s, r_pi, r_pipi])
        self.off = self.off + off_ext

    def checkRedudant(self,):
        terms = {}
        for i in self.angle:
            ang_term = "%02d_%02d_%02d"%(int(i[0]), int(i[1]), int(i[2]))
            if ang_term in terms.keys():
                print("Now:", i)
                print("Pre:", terms[ang_term])
            else:
                terms[ang_term] = i
    
    def checkout(self, atoms):
        b = Ffield()
        b.clearup()
        index = []
        
        b.lg = self.lg
        b.elements = atoms

        for i in self.atom:
            if i[0] in atoms:
                b.atom.append(i)

        for i in atoms:
            index.append(self.elements.index(i) + 1)
        
        for i in self.bond:
            i[0] = int(i[0])
            i[1] = int(i[1])
            if (i[0] in index) and  (i[1] in index):
                i[0] = index.index(i[0]) + 1
                i[1] = index.index(i[1]) + 1
                b.bond.append(i)

        for i in self.off:
            i[0] = int(i[0])
            i[1] = int(i[1])
            if (i[0] in index) and  (i[1] in index):
                i[0] = index.index(i[0]) + 1
                i[1] = index.index(i[1]) + 1
                b.off.append(i)
        
        for i in self.angle:
            i[0] = int(i[0])
            i[1] = int(i[1])
            i[2] = int(i[2])
            if (i[0] in index) and  (i[1] in index) and (i[2] in index):
                i[0] = index.index(i[0]) + 1
                i[1] = index.index(i[1]) + 1
                i[2] = index.index(i[2]) + 1
                b.angle.append(i)

        for i in self.torsion:
            i[0] = int(i[0])
            i[1] = int(i[1])
            i[2] = int(i[2])
            i[3] = int(i[3])
            if (i[0] in index) and  (i[1] in index) and (i[2] in index) and (i[3] in index):
                i[0] = index.index(i[0]) + 1
                i[1] = index.index(i[1]) + 1
                i[2] = index.index(i[2]) + 1
                i[3] = index.index(i[3]) + 1
                b.torsion.append(i)

        for i in self.hbond:
            i[0] = int(i[0])
            i[1] = int(i[1])
            i[2] = int(i[2])
            if (i[0] in index) and  (i[1] in index) and (i[2] in index):
                i[0] = index.index(i[0]) + 1
                i[1] = index.index(i[1]) + 1
                i[2] = index.index(i[2]) + 1
                b.hbond.append(i)
        return b

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
        # equation 3
        eq3 = []
        for i in range(len(self.atom)):
            val = int(float(self.atom[i][2]))
            val_boc = int(float(self.atom[i][28]))
            eq3.append([i, val, val_boc])

        # equation 4
        eq4 = []
        for i in range(len(self.atom)):
            for j in range(i, len(self.atom) - 1):
                p_boc1 = float(self.gl[0])
                p_boc2 = float(self.gl[1])
                p_o_131 = math.sqrt(float(self.atom[i][20])*float(self.atom[j][20]))
                p_o_132 = math.sqrt(float(self.atom[i][21])*float(self.atom[j][21]))
                p_o_133 = math.sqrt(float(self.atom[i][22])*float(self.atom[j][22]))
                eq4.append([i, j, p_boc1, p_boc2, p_o_132, p_o_131, p_o_133])
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

        # equation 7
        eq7 = []
        for i in range(len(self.atom)):
            val = int(float(self.atom[i][2]))
            val_e = int(float(self.atom[i][8]))
            nlp_opt = int(0.5 * (val_e - val))
            plp1 = float(self.gl[15])
            plp2 = float(self.atom[i][18])
            eq7.append([i, val_e, plp1, nlp_opt, plp2])
        
        # equation 11
        # equation 12
        # equation 13
        # equation 14
        # equation 15
        # equation 16
        # equation 17
        # equation 18

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

        # equation 24

        # equation 26
        eq26 = []
        for i in range(len(self.atom)):
            for j in range(i, len(self.atom) - 1):
                d = math.sqrt(float(self.atom[i][31])*float(self.atom[j][31]))
                #print(self.atom[i][0], self.atom[j][0], self.atom[i][31], self.atom[j][31], d
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
        
        o.write("# equation 3: coarse bond order\n") 
        o.write("%4s%8s%8s\n"%("a1", "val", "val_boc"))
        for i in range(len(eq3)):
            counter = 0 
            for j in range(len(eq3[i])):
                if counter < 1:
                    o.write("%4s"%self.elements[eq3[i][j]])
                else:
                    o.write("%8d"%eq3[i][j])
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

        o.write("# equation 7: lone pair energy\n") 
        o.write("%4s%8s%8s%8s%8s\n"%("a1", "val_e", "p_lp1", "nlp_opt", "p_lp2"))
        for i in range(len(eq7)):
            counter = 0 
            for j in range(len(eq7[i])):
                if counter < 1:
                    o.write("%4s"%self.elements[eq7[i][j]])
                else:
                    o.write("%9.4f"%eq7[i][j])
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

    def toParams(self,):
        """
        Generate params file for training
        """
        scale = 0.2
        o = open("params_total", "w")
        #@note: need to finish
              # 1        2         3         4      5
        GL = ["p_boc1", "p_boc2", "p_coa2", "Null", "Null", 
              # 6        7         8         9      10
              "Null", "p_ovun6", "Null", "Null", "Null", 
              # 11       12        13        14     15
              "Null", "Null", "Null", "Null", "Null", 
              # 16       17        18        19     20
              "p_lp1", "Null", "Null", "Null", "Null", 
              # 21       22        23        24     25
              "Null", "Null", "Null", "p_tor2", "p_tor2", 
              # 26       27        28        29     30
              "p_tor2", "Null", "p_cot2", "p_vdW1", "p_coa4", 
              # 31       32        33        34     35
              "Null", "p_ovun4", "p_ovun3", "Null", "Null", 
              # 36       37        38        39     40
              "Null", "Null", "Null", "p_coa3", "Null",]
              # 1        2         3         4      5
        GLE = ["4c", "4d", "15", "Null", "Null", 
              # 6        7         8         9      10
              "p_boc1", "12", "p_coa2", "Null", "Null", 
              # 11       12        13        14     15
              "p_boc1", "p_boc2", "p_coa2", "Null", "Null", 
              # 16       17        18        19     20
              "8", "p_boc2", "p_coa2", "Null", "Null", 
              # 21       22        23        24     25
              "p_boc1", "p_boc2", "p_coa2", "16a", "16b", 
              # 26       27        28        29     30
              "16b", "p_boc2", "17b", "23b", "15", 
              # 31       32        33        34     35
              "p_boc1", "11b", "11b", "Null", "Null", 
              # 36       37        38        39     40
              "p_boc1", "p_boc2", "p_coa2", "15", "Null",]
        n = 1
        for i in self.gl:
            val = i
            if val > 0:
                start = val * (1-scale)
                end = val * (1+scale)
            else:
                end = val * (1-scale)
                start = val * (1+scale)
            interval = abs(end -start) /20.0
            o.write("%4d%6d%12.4f%12.4f%12.4f%12.4f%12.4f"%(1, n, interval, start, end, val, val))
            o.write(" ! %s\n"%("gl"))
            n += 1

        ATOM = ["ro(sigma)", "Val", "mass", "Rvdw", "Dij", "gamma qeq", "ro(pi)", "Val(e)",
                "alfa", "gamma(w)", "Val(angle)", "p(ovun5)", "Null", "ChiEEM", "etaEEM", "Null",
                "ro(pipi)", "p(lp2)", "Heat", "p(boc4)", "p(boc3)", "p(boc5)", "Null", "Null",
                "p(ovun2)", "p(val3)", "Null", "Val(boc)", "p(val5)", "r_inner", "D_inner", "alpha_inner",
                "D_lg", "r_lg", ]
        ATOME = ["2", "3a", "Null", "23a", "23a", "24", "2", "7",
                 "23a", "23b", "13e", "12", "Null", "qeq", "qeq", "Null",
                 "2", "10", "Null", "4e", "4e", "4e", "Null", "Null",
                 "11a & 12", "13b", "Null", "3b", "13c", "25", "25", "25",
                 "26", "26",]
        n1 = 1
        for i in self.atom:
            n2 = 1
            for j in i[1:]:
                a1 = i[0]
                val = float(j)
                if val > 0:
                    start = val * (1-scale)
                    end = val * (1+scale)
                else:
                    end = val * (1-scale)
                    start = val * (1+scale)
                interval = abs(end -start) /20.0
                if n2 in [2, 3, 8, 11, 13, 16, 19, 23, 24, 27]:
                    pass
                else:
                    o.write("%4d%6d%6d%12.4f%12.4f%12.4f%12.4f"%(2, n1, n2, interval, start, end, val))
                    o.write(" ! %s %4s %s in %s\n"%("at", "@"+a1, ATOM[n2-1], ATOME[n2-1]))
                n2 += 1
            n1 += 1

        BOND = ["De(sigma)", "De(pi)", "De(pipi)", "p(be1)", "p(bo5)", "13corr", "p(bo6)", "p(ovun1)",
                "p(be2)", "p(bo3)", "p(bo4)", "Null", "p(bo1)", "p(bo2)", "ovc(3bond)", "Null"]
        BONDE = ["6", "6", "6", "6", "2", "Null", "2", "11a",  
                 "6", "2", "2", "Null", "2", "2", "Null", "Null", ]
        n1 = 1
        for i in self.bond:
            n2 = 1
            for j in i[2:]:
                a1 = self.elements[int(i[0]) -1]
                a2 = self.elements[int(i[1]) -1]
                val = float(j)
                start = val * (1-scale)
                end = val * (1+scale)
                interval = abs(end -start) /20.0
                if n2 in [6, 12, 15, 16]:
                    pass
                else:
                    o.write("%4d%6d%6d%12.4f%12.4f%12.4f%12.4f"%(3, n1, n2, interval, start, end, val))
                    o.write(" ! %s %4s %4s %s in %s\n"%("bo", "@"+a1, "@"+a2, BOND[n2-1], BONDE[n2-1]))
                n2 += 1
            n1 += 1

        OFF = ["Dij", "RvdW", "alfa", "ro(sigma)", "ro(pi)", "ro(pipi)", "D_lg"]
        OFFE = ["23a", "23a", "23a", "2", "2", "2", "26"]
        n1 = 1
        for i in self.off:
            n2 = 1
            for j in i[2:]:
                a1 = self.elements[int(i[0]) -1]
                a2 = self.elements[int(i[1]) -1]
                val = float(j)
                start = val * (1-scale)
                end = val * (1+scale)
                interval = abs(end -start) /20.0
                o.write("%4d%6d%6d%12.4f%12.4f%12.4f%12.4f"%(4, n1, n2, interval, start, end, val))
                o.write(" ! %s %4s %4s %s in %s\n"%("of", "@"+a1, "@"+a2, OFF[n2-1], OFFE[n2-1]))
                n2 += 1
            n1 += 1

        ANG = ["Theta", "p_val1", "p_val2", "p_coa1", "p_val7", "p_pen1", "p_val4"]
        ANGE = ["13g", "13a", "13a", "15", "13c", "14a", "13b"]
        n1 = 1
        for i in self.angle:
            n2 = 1
            for j in i[3:]:
                a1 = self.elements[int(i[0]) -1]
                a2 = self.elements[int(i[1]) -1]
                a3 = self.elements[int(i[2]) -1]
                val = float(j)
                start = val * (1-scale)
                end = val * (1+scale)
                interval = abs(end -start) /20.0
                o.write("%4d%6d%6d%12.4f%12.4f%12.4f%12.4f"%(5, n1, n2, interval, start, end, val))
                o.write(" ! %s %4s %4s %4s %s in %s\n"%("ang", "@"+a1, "@"+a2, "@"+a3, ANG[n2-1], ANGE[n2-1]))
                n2 += 1
            n1 += 1

        TOR = ["V1", "V2", "V3", "p_tor1", "p_cot1", "Null", "Null"]
        TORE = ["16a", "16a", "16a", "16a", "17a", "", ""]
        n1 = 1
        for i in self.torsion:
            n2 = 1
            for j in i[4:]:
                a1 = self.elements[int(i[0]) -1]
                a2 = self.elements[int(i[1]) -1]
                a3 = self.elements[int(i[2]) -1]
                a4 = self.elements[int(i[3]) -1]
                val = float(j)
                start = val * (1-scale)
                end = val * (1+scale)
                interval = abs(end -start) /20.0
                if n2 in [5, 6]:
                    pass
                else:
                    o.write("%4d%6d%6d%12.4f%12.4f%12.4f%12.4f"%(6, n1, n2, interval, start, end, val))
                    o.write(" ! %s %4s %4s %4s %4s %s in %s\n"%("tor", "@"+a1, "@"+a2, "@"+a3, "@"+a4, TOR[n2-1], TORE[n2-1]))
                n2 += 1
            n1 += 1

        HBO = ["r_hb", "p_hb1", "p_hb2", "p_hb3"]
        HBOE = ["18", "18", "18", "18"]
        n1 = 1
        for i in self.hbond:
            n2 = 1
            for j in i[3:]:
                a1 = self.elements[int(i[0]) -1]
                a2 = self.elements[int(i[1]) -1]
                a3 = self.elements[int(i[2]) -1]
                val = float(j)
                if val > 0:
                    start = val * (1-scale)
                    end = val * (1+scale)
                else:
                    end = val * (1-scale)
                    start = val * (1+scale)
                interval = abs(end -start) /20.0
                o.write("%4d%6d%6d%12.4f%12.4f%12.4f%12.4f"%(5, n1, n2, interval, start, end, val))
                o.write(" ! %s %4s %4s %4s %s in %s\n"%("hbo", "@"+a1, "@"+a2, "@"+a3, HBO[n2-1], HBOE[n2-1]))
                n2 += 1
            n1 += 1

        o.close()
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", default="ffield", nargs="?", help="force field file name")
    parser.add_argument("-D", action="store_true", help="Debug the code")
    parser.add_argument("-trans", action="store_true", help="transform ffield to more readable format")
    parser.add_argument("-params", action="store_true", help="generate params for training")
    parser.add_argument("-complete", action="store_true", help="complete the off table")
    parser.add_argument("-type", nargs=1, type=int, help="Force field type: 0 for vdw; 1 for lg_inner wall")
    parser.add_argument("-checkout", nargs='+', help="check out the force field")
    parser.add_argument("-check", action="store_true", help="check the force field")
    args = parser.parse_args()
    #print(b.getBondDist(3,2)
    
    fname = args.fname

    assert os.path.exists(fname)

    if args.D:
        DEBUG = 1

    if args.type:
        ntype = args.type[0]
    else:
        ntype = 0
        print("Note: Using default force field type")

    ff = Ffield(fname, ntype)
    if args.trans:
        ff.toEquation()
    if args.params:
        ff.toParams()
    if args.complete:
        ff.completeOff()
        toFfield(ff)
    if args.checkout:
        atoms = args.checkout
        ff2 = ff.checkout(atoms)
        toFfield(ff2)
    if args.check:
        ff.checkRedudant()

def test1():
    fname = "ffield"
    ntype = 0
    ff = Ffield(fname, ntype)
    ff2 = ff.checkout()
    toFfield(ff2)

def test():
    fname = "ffield"
    ntype = 0
    ff = Ffield(fname, ntype)
    ff.checkRedudant()
    
if __name__ == "__main__":
    #test()
    main()
