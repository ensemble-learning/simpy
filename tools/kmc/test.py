T = 300.0
class Rxn():
    def __init__(self,):
        self.rxns = []
        self.species = []

def read_rxn(rxn):
    delta_e, delta_s = [], []
    delta_g = []
    delta_e_act, delta_s_act = [], []
    delta_g_act = []
    species = []
    is_sur = []
    f = open('rxn.csv', 'r')
    for i in f:
        if not i.strip().startswith('#'):
            tokens = i.strip().split(',')
            delta_e.append(float(tokens[1].strip()))
            delta_s.append(float(tokens[2].strip()))
            delta_e_act.append(float(tokens[3].strip()))
            delta_s_act.append(float(tokens[4].strip()))
            r = tokens[0].split('<-->')[0]
            p = tokens[0].split('<-->')[1]
            reactants = [ii.strip() for ii in r.split('+')]
            for ii in reactants:
                if not ii in species:
                    species.append(ii)
            products = [ii.strip() for ii in p.split('+')]
            for ii in products:
                if not ii in species:
                    species.append(ii)
    f.close()
    print(species)
    for e, s in zip(delta_e, delta_s):
        delta_g.append(e-T*s)
    for e, s in zip(delta_e_act, delta_s_act):
        delta_g_act.append(e-T*s)


rxn = Rxn()
read_rxn(rxn)


