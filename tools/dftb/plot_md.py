"""Plot the pressure, potential energy and temperature of the simulation
"""

import matplotlib.pyplot as plt
import numpy as np

def plt_all():
    step = []
    press = []
    pot = []
    total = []
    tempt = []
    
    f = open("md.out", "r")
    
    for i in f:
        tokens = i.strip().split(":")
        if "MD step" in i:
            val = tokens[1].split()[0]
            step.append(int(val))
        if "Pressure" in i:
            val = tokens[1].split()[2]
            press.append(float(val))
        if "Potential Energy" in i:
            val = tokens[1].split()[2]
            pot.append(float(val))
        if "Total MD Energy" in i:
            val = tokens[1].split()[2]
            total.append(float(val))
        if "Temperature" in i:
            val = tokens[1].split()[2]
            tempt.append(float(val))
         
    f.close()
    
    step = np.array(step)
    press = np.array(press)
    pot = np.array(pot)
    tempt = np.array(tempt)
    
    fig = plt.figure()
    
    ax = fig.add_subplot(4,1,1)
    ax.set_title("Pressure")
    ax.set_ylabel("Pa")
    if len(step) > len(press):
        x = step[:len(press)]
    ax.plot(x, press)
    np.savetxt("pressure", np.array([x, press]))
    
    ax = fig.add_subplot(4,1,2)
    ax.set_title("Potential Energy")
    ax.set_ylabel("ev")
    if len(step) > len(pot):
        x = step[:len(pot)]
    ax.plot(x, pot)
    np.savetxt("potential", np.array([x, pot]))
    
    ax = fig.add_subplot(4,1,3)
    ax.set_title("Temperature")
    ax.set_ylabel("K")
    if len(step) > len(tempt):
        x = step[:len(tempt)]
    ax.plot(x, tempt)
    np.savetxt("temperature", np.array([x, tempt]))
    
    ax = fig.add_subplot(4,1,4)
    ax.set_title("Total Energy")
    ax.set_ylabel("ev")
    if len(step) > len(total):
        x = step[:len(total)]
    ax.plot(x, total)
    np.savetxt("total", np.array([x, total]))
    
    plt.subplots_adjust(hspace=0.5)
    plt.savefig("md.png")
    plt.show()
    
def plt_no_pressure():
    step = []
    pot = []
    total = []
    tempt = []
    
    f = open("md.out", "r")
    
    for i in f:
        tokens = i.strip().split(":")
        if "MD step" in i:
            val = tokens[1].split()[0]
            step.append(int(val))
        if "Potential Energy" in i:
            val = tokens[1].split()[2]
            pot.append(float(val))
        if "Total MD Energy" in i:
            val = tokens[1].split()[2]
            total.append(float(val))
        if "Temperature" in i:
            val = tokens[1].split()[2]
            tempt.append(float(val))
         
    f.close()
    
    step = np.array(step)
    pot = np.array(pot)
    tempt = np.array(tempt)
    
    fig = plt.figure()
    
    ax = fig.add_subplot(3,1,1)
    ax.set_title("Potential Energy")
    ax.set_ylabel("ev")
    if len(step) > len(pot):
        x = step[:len(pot)]
    else:
        x = step
    ax.plot(x, pot)
    output = np.array([x, pot])
    output = output.transpose()
    np.savetxt("potential", output)
    
    ax = fig.add_subplot(3,1,2)
    ax.set_title("Temperature")
    ax.set_ylabel("K")
    if len(step) > len(tempt):
        x = step[:len(tempt)]
    else:
        x = step
    ax.plot(x, tempt)
    output = np.array([x, tempt])
    output = output.transpose()
    np.savetxt("temperature", output)
    
    ax = fig.add_subplot(3,1,3)
    ax.set_title("Total Energy")
    ax.set_ylabel("ev")
    if len(step) > len(total):
        x = step[:len(total)]
    else:
        x = step
    ax.plot(x, total)
    output = np.array([x, total])
    output = output.transpose()
    np.savetxt("total", output)
    
    plt.subplots_adjust(hspace=0.5)
    plt.savefig("md.png")
    plt.show()
    

if __name__ == "__main__":
    plt_no_pressure()
