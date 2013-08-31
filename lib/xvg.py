import numpy as np
import matplotlib.pyplot as plt
import re

class Xvg():
    def __init__(self, filename):
        self.data = np.array([])
        self.notes = []
        self.description = []
        self.title = ""
        self.xlabel = ""
        self.ylabel = ""
        self.xscale = 1.0
        self.yscale = 1.0
        self.ny = 0
        self.legend = []
        self.read(filename)
    def read(self, filename):
        data = []
        f = open(filename, "r")
        for i in f:
            if i.strip().startswith("#"):
                self.notes.append(i)
            elif i.strip().startswith("@"):
                self.description.append(i)
            else:
                tokens = i.strip().split()
                data.append([float(j) for j in tokens])
        data = np.array(data)
        data = data.transpose()
        self.data = data
    def parser(self, ):
        p = re.compile('.*"(.*)"')       
        for i in self.description:
            if "title" in i:
                m = p.match(i)
                if m:
                    self.title = m.group(1)
            elif "xaxis" in i:
                m = p.match(i)
                if m:
                    self.xlabel = m.group(1)
            elif "yaxis" in i:
                m = p.match(i)
                if m:
                    self.ylabel = m.group(1)
            elif "@ s" in i:
                m = p.match(i)
                if m:
                    self.legend.append(m.group(1))
                self.ny += 1
                    
    def plot(self, ):
        plt.plot(self.data[0], self.data[1])
        plt.title(self.title, size="x-large")
        plt.xlabel(self.xlabel, size="x-large")
        plt.ylabel(self.ylabel, size="x-large")
        if 1:
            aver, std = self.statistic(1)
            plt.plot(self.data[0], [aver]*len(self.data[0]), 
            ls="--", lw=4, color="black")
        plt.show()
    
    def statistic(self, n):
        start = int(len(self.data[n]) * 0.2)
        aver = np.average(self.data[n][start:])
        std = np.std(self.data[n][start:])
        print aver*2, std*2
        return aver, std
        
xvg = Xvg("area.xvg")
xvg.parser()
xvg.statistic(1)
xvg.plot()
        
