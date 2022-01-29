class ML_AB():
    def __init__(self,):
        self.n_config = 0
        self.n_atom_typte = 0
        self.n_atoms = 0
        self.pe_ref = []
        self.masses = []
        self.n_basis_sets = []
        self.head = []
        self.confs = []

    def read(self, fname='ML_AB'):
        f = open('ML_ABN', 'r')
        lines = f.readlines()
        f.close()

        blocks = []
        blocks.append([])
        nb = 0
        for i in lines:
            tokens = i.strip()
            if tokens.startswith('Configuration num.'):
                blocks.append([])
                nb+=1
            blocks[nb].append(i)
        self.head = blocks[0]
        self.confs = blocks[1:]
        
        self.n_config = int(lines[4].strip())
        self.n_atom_typte = int(lines[8].strip())
        self.n_basis_sets = [int(ii) for ii in lines[32].strip().split()]

    def write(self,):
        o = open('ml-ab-head.dat', 'w')
        for i in self.head:
            o.write(i)
        o.close()
        o = open('ml-ab-confs.dat', 'w')
        for i in self.confs:
            for j in i:
                o.write(j)
        o.close()
        o = open('ml-ab-log.dat', 'w')
        o.write('%d\n'%self.n_config)
        o.write(' '.join(['%d'%ii for ii in self.n_basis_sets]))
        o.write('\n')
        o.close()


ab = ML_AB()
ab.read()
ab.write()
