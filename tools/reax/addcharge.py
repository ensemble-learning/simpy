"""
convert qm charge to reaxFF charge training set
"""

import os

def main():
    if os.path.exists("qm_charge.dat"):
        o = open("add.charge", "w")
        f = open("qm_charge.dat", "r")
        for i in f:
            tokens = i.strip().split()
            if len(tokens) > 4:
                mol = tokens[0]
                n = int(tokens[2])
                q = float(tokens[3])
                o.write("%-20s%-8.4f%5d%12.4f\n"%(mol, 0.1, n, q))
        o.close()

if __name__ == "__main__":
    main()
