"""
Manipulate the csv file
"""

import csv
import numpy as np
import matplotlib.pyplot as plt

class Csv():
    def __init__(self, fname):
        self.n = 0
        self.data = []
        self.name = fname
        self.read()
    def read(self,):
        data = []

        f = open("H2.csv", "r")

        reader = csv.reader(f)

        for row in reader:
            data.append([float(j) for j in row])

        data = np.array(data)
        data = data.transpose()
        
        self.data = data

def test():
    a = Csv("H2.csv")
    plt.plot(a.data[0], a.data[1])
    plt.show()

if __name__ == "__main__":
    test()
    
