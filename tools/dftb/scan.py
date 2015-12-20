import numpy as np
import matplotlib.pyplot as plt

f = open("geo_start.gen", "r")
lines = f.readlines()
tokens = lines[0].strip().split()
natoms = int(tokens[0])


cells = []

for i in range(natoms+3, natoms+6):
    tokens = lines[i].strip().split()
    cells.append([float(j) for j in tokens])


data = np.linspace(0.8, 1.2, 11)
data = np.power(data, 1/3.0)
for i in range(len(data)):
    print data
    
    

