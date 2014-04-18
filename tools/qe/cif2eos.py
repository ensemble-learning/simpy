#!/home/liulc/src/epd-7.1/epd-7.1-2-rh5-x86/bin/python
# -*- coding: utf-8 -*-

import os, re, sys
from sys import argv
from math import cos, sin, radians, pow

mass = {'H':1.00794, 'He':4.002602,'Li':6.941,  'Be':9.012182,\
        'B':10.811, 'C':12.0107, 'N':14.00674,'O':15.9994, \
        'F':18.9984032, 'Ne':20.1797,'Na':22.98977,'Mg':24.3050,\
        'Al':26.981538, 'Si':28.0855,'P':30.973761,'S':32.066, \
        'Cl':35.4527, 'Ar':39.948,'K':39.0983,'Ca':40.078, \
		'Pd':106.42, 'Mo': 95.94, 'Nb':92.91, 'Te':127.6, \
		'V':50.94}

ang2bohr = 1.889725989

TEMPLATE_EOS = """\
  &CONTROL
                        title = $JOBNAME$ ,
                  calculation = 'relax' ,
                 restart_mode = 'from_scratch' ,
                   pseudo_dir = '$PSEUDO_PATH$' ,
                       prefix = '$JOBNAME$' ,
  /
  &SYSTEM
                        ibrav = 14,
                    celldm(1) = $A$,  
                    celldm(2) = $B/A$,  
                    celldm(3) = $C/A$,
                    celldm(4) = $COS_ALPHA$, 
                    celldm(5) = $COS_BETA$,
                    celldm(6) = $COS_GAMMA$,
                          nat = $NATOMS$,
                         ntyp = $NTYPE$,
                      ecutwfc = 30 ,
                      ecutrho = 240 , 
#                        nosym = .true. ,
                  occupations = 'smearing' ,
                      degauss = 0.01 ,
                     smearing = 'marzari-vanderbilt' ,
                        nspin = 1 ,
  /
  &ELECTRONS
             electron_maxstep = 500,
                  mixing_beta = 0.2 ,
                  mixing_mode = 'local-TF' ,
                  mixing_ndim = 25,
  /
 # &IONS
 #                ion_dynamics = 'bfgs' ,
 #           pot_extrapolation = 'second_order' ,
 #           wfc_extrapolation = 'second_order' ,
 # /
 # &CELL
 #               cell_dynamics = 'bfgs' ,
 #                 cell_dofree = 'xyz' ,
 # /
 ATOMIC_SPECIES
$ATOMINFO$
 ATOMIC_POSITIONS crystal 
$XYZ$
 K_POINTS automatic 
   2 2 2   0 0 0   
"""

class CIF_DATA:
	def __init__(self, filename):
		self.read_cif(filename)

	def read_cif(self, filename):
		cif_file = open(filename, 'r')
		xyz_tokens = []
		for i in cif_file:
			#Box information
			if re.search("length_a", i):    a    = ang2bohr * float(i.split()[1])
			if re.search("length_b", i):    b    = ang2bohr * float(i.split()[1])
			if re.search("length_c", i):    c    = ang2bohr * float(i.split()[1])
			if re.search("angle_alpha", i): alpha = float(i.split()[1])
			if re.search("angle_beta", i):  beta  = float(i.split()[1])
			if re.search("angle_gamma", i): gamma = float(i.split()[1])
			#relative xyz coord
			if re.search("Uiso", i): 
				xyz_tokens.append("%-3s%10s%10s%10s" %(i.split()[1], i.split()[2], i.split()[3], i.split()[4]))
			if re.search("Uani", i): 
				xyz_tokens.append("%-3s%10s%10s%10s" %(i.split()[1], i.split()[2], i.split()[3], i.split()[4]))
		cif_file.close()
		self.a = a
		self.b2a = str(b/a)
		self.c2a = str(c/a)
		self.cos_a = "%.8f" %cos(radians(alpha))
		self.cos_b = "%.8f" %cos(radians(beta))
		self.cos_c = "%.8f" %cos(radians(gamma))
		self.xyz = "\n".join(xyz_tokens)
		#check natoms, ntype
		self.natoms = str(len(xyz_tokens))
		self.atom_dic = {}
		for i in xyz_tokens:
			self.atom_dic[i.split()[0]] = 1
		self.ntype = str(self.atom_dic.__len__())
		atominfo_tokens = []
		for i in self.atom_dic.keys():
			atominfo_tokens.append("%-5s%10.5f%5s.pbe-rrkjus.UPF" %(i, mass[i], i))
		self.atominfo = "\n".join(atominfo_tokens)

def replace_temp(c, fout_name, jobname, ratio=1.0):
	fout = open(fout_name, 'w')
	text1 = re.sub("\$JOBNAME\$", jobname, TEMPLATE_EOS)
	text2 = re.sub("\$PSEUDO_PATH\$", '/home/liulc/pseudo/', text1)
	text3 = re.sub("\$A\$", str(c.a*pow(ratio, 1./3.)), text2) #ratio is the scale ratio
	text4 = re.sub("\$B/A\$", c.b2a, text3)
	text5 = re.sub("\$C/A\$", c.c2a, text4)
	text6 = re.sub("\$COS_ALPHA\$", c.cos_a, text5)
	text7 = re.sub("\$COS_BETA\$", c.cos_b, text6)
	text8 = re.sub("\$COS_GAMMA\$", c.cos_c, text7)
	text9 = re.sub("\$NATOMS\$", c.natoms, text8)
	text10 = re.sub("\$NTYPE\$", c.ntype, text9)
	text11 = re.sub("\$ATOMINFO\$", c.atominfo, text10)
	text12 = re.sub("\$XYZ\$", c.xyz, text11)
	text13 = re.sub("Si.pbe-rrkjus.UPF", "Si.pbe-rrkj.UPF", text12)
	fout.write(text13)
	fout.close()

if __name__ == "__main__":
	c = CIF_DATA(argv[1])
	for i in [j*0.01 for j in range(80, 145, 5)]:  #the eos range according to volume
		replace_temp(c, argv[1].split('.')[0]+"_"+str(i)+".in", argv[1].split('.')[0]+"_"+str(i), i)
		cmd = "exQespresso "+argv[1].split('.')[0]+"_"+str(i)
		#os.system(cmd)
	
