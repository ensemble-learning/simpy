import time
import os

def msd2pdb(file_name):
    f = open(file_name + ".msd" ,'r')
    time_now = time.strftime('%Y-%m-%d',time.localtime(time.time()))
    o = open(file_name + "_msd.pdb", 'w')
    label = ["ATOM"] 

    o.write("""HEADER %s 
REMARK    MADE BY DIRECT FORCE FIELD.
REMARK %s 
"""%(file_name, time_now))
# msd file coordinate
# 1     C0    6   CCH3        c_4    0.000000  0.000000    0.000000    0.765427   1    UNK  0
    for i in f:
        if len(i.split()) == 12:
            a = i.split()
            atom_name = a[1][0] + a[0]
            res_no = int(a[9])
            res_name = a[10]
            o.write("%-6s%5d  %-3s %3s  %4d  %8.3f%8.3f%8.3f\n"%(label[0],int(a[0]), atom_name, res_name, res_no, float(a[6]), float(a[7]), float(a[8])))
    o.close()
    f.close()

if __name__ == "__main__":
    for i in os.listdir("."):
        if i.endswith(".msd"):
            msd2pdb(i.split('.')[0])
    
