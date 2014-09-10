"""
0: step
1: total energy   
2: poten. energy       
3: amd E   
4: delta. energy    
5: force. scale     
6: kin. energy      
7: temp.     
8: target       
9: volume       
10: press.       
11:target
"""
import numpy as np
import os

data = np.loadtxt("water.out", skiprows=1)
data = data.transpose()
np.savetxt("tempt.txt", data[7].transpose(), fmt="%.4f")
