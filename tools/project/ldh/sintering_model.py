""" A simple model to test the idea of surface reducation
when two (or more) sphere merge into one large sphere
"""

import math
def cal(n = 12):
    r1 = 4.0
    r2 = math.pow(0.2, 1/3.0) * r1
    r3 = math.pow(n, 1/3.0) * r1
    r4 = r3 - (r1 - r2)
    ratio1 = 1 - math.pow(r4/r3, 3.0)
    
    v = n * 4.0/3 * math.pi * r1**3
    a1 = math.pow(v, 1/3.0)
    a2 = a1 - (r1 - r2)
    ratio2 = 1 - math.pow(a2/a1, 3.0)
    print "%4d, %8.2f, %8.2f"%(n, ratio1, ratio2)

for i in range(1, 100):
    cal(i)