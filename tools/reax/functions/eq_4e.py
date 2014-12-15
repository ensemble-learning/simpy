import numpy as np
import matplotlib.pyplot as plt

def f4(delta_boc, p_boc3, p_boc4, p_boc5):
    bo = np.linspace(0, 1.0, 41)

    e1 = -p_boc3*(p_boc4*bo*bo - delta_boc) + p_boc5
    scale = 1/(1 + np.exp(e1))
    
    return bo, scale


p_boc3 = 2.8410
p_boc4 = 4.104
p_boc5 = 0.0003

delta_boc =  0.0

fig = plt.figure()
ax = fig.add_subplot(221)
delta_boc = -0.2
x, y = f4(delta_boc, p_boc3, p_boc4, p_boc5)
ax.plot(x, y, label="decrease", lw=3)


delta_boc =  0.0
x, y = f4(delta_boc, p_boc3, p_boc4, p_boc5)
ax.plot(x, y, label="ref", lw=3)

delta_boc =  0.2
x, y = f4(delta_boc, p_boc3, p_boc4, p_boc5)
ax.plot(x, y, label="increase", lw=3)
ax.legend()

ax = fig.add_subplot(222)
p_boc3 = 2.841 * 0.8
x, y = f4(delta_boc, p_boc3, p_boc4, p_boc5)
ax.plot(x, y, label="decrease", lw=3)

p_boc3 = 2.841 * 1.0
x, y = f4(delta_boc, p_boc3, p_boc4, p_boc5)
ax.plot(x, y, label="ref", lw=3)

p_boc3 = 2.841 * 1.2
x, y = f4(delta_boc, p_boc3, p_boc4, p_boc5)
ax.plot(x, y, label="increase", lw=3)
ax.legend()

ax = fig.add_subplot(223)
p_boc4 = 4.104 * 0.8
x, y = f4(delta_boc, p_boc3, p_boc4, p_boc5)
ax.plot(x, y, label="decrease", lw=3)

p_boc4 = 4.104 * 1.0
x, y = f4(delta_boc, p_boc3, p_boc4, p_boc5)
ax.plot(x, y, label="ref", lw=3)

p_boc4 = 4.104 * 28.2
x, y = f4(delta_boc, p_boc3, p_boc4, p_boc5)
ax.plot(x, y, label="increase", lw=3)
ax.legend()

ax = fig.add_subplot(224)
p_boc5 = 1 * 0.8
x, y = f4(delta_boc, p_boc3, p_boc4, p_boc5)
ax.plot(x, y, label="decrease", lw=3)

p_boc5 = 1 * 1.0
x, y = f4(delta_boc, p_boc3, p_boc4, p_boc5)
ax.plot(x, y, label="ref", lw=3)

p_boc5 = 1 * 2.2
x, y = f4(delta_boc, p_boc3, p_boc4, p_boc5)
ax.plot(x, y, label="increase", lw=3)
ax.legend()

plt.show()
