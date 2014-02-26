import numpy as np
import matplotlib.pyplot as plt

p_bo1 = -0.10
p_bo2 = 8.0
r_o = 2.1

r = np.linspace(1.0, 6.0, 41)
bo = np.exp(p_bo1*np.power(r/r_o, p_bo2))

plt.plot(r, bo)
plt.show()
