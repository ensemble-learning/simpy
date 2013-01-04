""" Fortran reactive force field Fort file
"""

class Fort73():
    """ data in fort.73 (energy components)
    """
    def __init__(self, filename='fort.73'):
        self.iter = []
        self.ebond = []
        self.eatom = []
        self.elp = []
        self.emol = []
        self.eval = []
        self.ecoa = []
        self.ehbo = []
        self.etors = []
        self.econj = []
        self.evdw = []
        self.ecol = []
        self.echarge = []
        self.efield = []
        self.read(filename)

    def read(self, filename):
        """ read data from fort.73
        """
        f = open(filename, 'r')
        counter = 1
        for i in f:
            if counter % 3 == 0:
                tokens = i.strip().split()
                self.iter.append(int(tokens[0]))
                self.ebond.append(float(tokens[1]))
                self.eatom.append(float(tokens[2]))
                self.elp.append(float(tokens[3]))
                self.emol.append(float(tokens[4]))
                self.eval.append(float(tokens[5]))
                self.ecoa.append(float(tokens[6]))
                self.ehbo.append(float(tokens[7]))
                self.etors.append(float(tokens[8]))
                self.econj.append(float(tokens[9]))
                self.evdw.append(float(tokens[10]))
                self.ecol.append(float(tokens[11]))
                self.echarge.append(float(tokens[12]))
                self.efield.append(float(tokens[13]))
            counter += 1

        f.close()

if __name__ == "__main__":
    print __doc__
