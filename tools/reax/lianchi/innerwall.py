#!/home/liulc/src/epd-7.1/epd-7.1-2-rh5-x86/bin/python
#-*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import math

# funtion of innerwall
# d is well depth, ro is equilibrium distance, alpha is alpha....
def innerwall(d,ro,alpha,r):
	return d*math.exp(alpha*(1-r/ro))

r = np.arange(0.4,5,0.01)
d = 0.1000
ro = 1.6737
alpha = 12.0000
y,y2 = [],[] 

for i in r:
	y.append(innerwall(d,ro,alpha,i))
for i in r:
	y2.append(innerwall(d,2.0,alpha,i))

plt.plot(r,y,'r-', label="r=1.67")
plt.plot(r,y2,'b-', label="r=2.0")
plt.legend()
plt.xlim(1,3.0)
plt.ylim(-1,130.0)
plt.show()

