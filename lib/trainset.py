"""
"""
import os
import numpy as np

class Trainset():
    def __init__(self,):
        self.energy = []
        self.energy_lines = []
        self.parse_fort99()

    def parse_fort99(self,):
        assert os.path.exists("fort.99")
        f = open("fort.99", 'r')
        reax = []
        qm = []

        for i in f:
            if i.strip().startswith("Energy"):
                self.energy_lines.append(i)
                reax.append(float(i[61:72]))
                qm.append(float(i[73:84]))
                break

        for i in f:
            if len(i.strip()) > 84:
                self.energy_lines.append(i)
                reax.append(float(i[61:72]))
                qm.append(float(i[73:84]))

        f.close()
        reax = np.array(reax)
        qm = np.array(qm)
        self.energy = [reax, qm]

def test():
    a = Trainset()
    
if __name__ == "__main__":
    test()
