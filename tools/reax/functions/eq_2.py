import numpy as np
import matplotlib.pyplot as plt

def eq2(r):
    p_bo1 = -0.10
    p_bo2 = 8.0
    p_bo3 = 1.0
    r_o = 2.1

    C12 = p_bo1 * np.power(r/r_o, p_bo2)
    bo = p_bo3 * np.exp(C12)
    bo_dev = p_bo3 * p_bo1/r_o * bo * p_bo2 * np.power(r/r_o, p_bo2 -1 )
    
    return bo, bo_dev


def eq6(bo, bo_dev):
    De = 20.0
    pbe1 = 1.0
    pbe2 = 1.0
    pbe3 = 1.0
    exp1 = np.exp(pbe1*(pbe3-np.power(bo, pbe2)))
    e_bond = -De*bo*exp1
    
    C1 = -De * exp1
    C2 = -e_bond * pbe1 * pbe2 * np.power(pbe3-bo, pbe2 -1)
    e_bond_dev = (C1 + C2)*bo_dev
    return e_bond, e_bond_dev

r = np.linspace(1.0, 6.0, 41)
bo, bo_dev = eq2(r)
e_bond, e_bond_dev = eq6(bo, bo_dev)

"""
plt.plot(r, bo_dev, lw=2)
plt.plot(r, bo, lw=2)
"""
plt.plot(r, e_bond, lw=2)
plt.plot(r, e_bond_dev, lw=2)
plt.show()

for i in range(len(r)):
    print r[i], e_bond[i], e_bond_dev[i]
