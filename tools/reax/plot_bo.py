"""
plot the bond order curve
Usage: python plot_bo.py atom1 atom2 r0
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_bo(atom1, atom2, r0):

    params = []
    
    f = open("output.ff", "r")
    
    for i in f:
        if "bond order" in i:
            break
    
    for i in f:
        if i.strip().startswith("#"):
            break
        else:
            tokens = i.strip().split()
            params.append(tokens)
    f.close()
    
    for i in params:
        a1 = i[0]
        a2 = i[1]
        if (atom1 == a1 and atom2 == a2) or (atom1 == a2 and atom2 == a1):
            ros = float(i[2])
            pbo1 = float(i[3])
            pbo2 = float(i[4])
            rop = float(i[5])
            pbo3 = float(i[6])
            pbo4 = float(i[7])
            ropp = float(i[8])
            pbo5 = float(i[9])
            pbo6 = float(i[10])
    
    x = np.linspace(0.3*ros, 4*ros, 100)
    ys = np.exp(pbo1*np.power(x/ros, pbo2))
    yt = ys
    plt.plot(x, ys)
    
    flag = 0
    if rop > 0.001:
        yp = np.exp(pbo3*np.power(x/rop, pbo4))
        yt += yp
        plt.plot(x, yp)
        flag += 1
    
    if ropp > 0.001:
        ypp = np.exp(pbo5*np.power(x/ropp, pbo6))
        yt += ypp
        plt.plot(x, ypp)
        flag += 1
    
    if flag:
        plt.plot(x, yt)
    
    
    plt.plot([r0, r0], [0, np.max(yt)], ls="--")
    bo1 = np.exp(pbo1*np.power(r0/ros, pbo2))
    bot = bo1
    if rop > 0.001:
        bo2 = np.exp(pbo3*np.power(r0/rop, pbo4))
        bot += bo2
    
    if ropp > 0.001:
        bo3 = np.exp(pbo5*np.power(r0/ropp, pbo6))
        bot += bo3
    
    plt.text(r0+0.1, 0.1, "r$_0$ = %.4f"%r0)
    
    if flag == 0:
        plt.text(r0+0.1, bot, "BO$^{\sigma}$ = %.4f"%bo1)
    elif flag == 1:
        plt.text(r0+0.1, bot+0.1, "BO$^{\sigma}$ = %.4f\nBO$^{\pi}$ = %.4f\nBO = %.4f"%(bo1, bo2, bot))
    elif flag == 2:
        plt.text(r0+0.1, bot+0.1, "BO$^{\sigma}$ = %.4f\nBO$^{\pi}$ = %.4f\nBO$^{\pi\pi}$ = %.4f\nBO = %.4f"%(bo1, bo2, bo3, bot))
    
    plt.xlabel("%s-%s distance ($\AA$)"%(atom1, atom2))
    plt.ylabel("Bond Order")
    plt.ylim([0, np.max(yt) + 0.2])
    
    n = 0
    for i in yt:
        n += 1
        if i < 0.0001:
            break
    plt.xlim([x[0], x[n]])
    
    # print the parameters
    
    if flag == 0:
        plt.text(r0*0.5 + x[n]*0.5, bot-0.3, "r$^{\sigma}$ = %.4f\np$_{bo1}$ = %.4f\np$_{bo2}$ = %.4f"%(ros, pbo1, pbo2))
    elif flag == 1:
        plt.text(r0*0.5 + x[n]*0.5, bot-0.5, "r$^{\sigma}$ = %.4f\np$_{bo1}$ = %.4f\np$_{bo2}$ = %.4f\n\
    r$^{\pi}$ = %.4f\np$_{bo3}$ = %.4f\np$_{bo4}$ = %.4f\n\
    "%(ros, pbo1, pbo2, rop, pbo3, pbo4))
    elif flag == 2:
        plt.text(r0*0.5 + x[n]*0.5, bot-0.7, "r$^{\sigma}$ = %.4f\np$_{bo1}$ = %.4f\np$_{bo2}$ = %.4f\n\
    r$^{\pi}$ = %.4f\np$_{bo3}$ = %.4f\np$_{bo4}$ = %.4f\nr$^{\pi\pi}$ = %.4f\np$_{bo5}$ = %.4f\np$_{bo6}$ = %.4f\n\
    "%(ros, pbo1, pbo2, rop, pbo3, pbo4, ropp, pbo5, pbo6))
            
    plt.show()

def main():
    if len(sys.argv) < 3:
        print __doc__
    else:
        atom1 = sys.argv[1]
        atom2 = sys.argv[2]
        r0 = float(sys.argv[3])
        plot_bo(atom1, atom2, r0)

if __name__ == "__main__":
    main()
