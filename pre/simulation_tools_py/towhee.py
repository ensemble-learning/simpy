"""transfer msd file to towhee input
@version: 1.0
@author: hawkweedcheng
@contact: chengtao@sjtu.edu.cn
"""
from mol import Atom, Bond, Molecule, System
import random

class Towhee():
    """Basic class for towhee files
    @todo: 
    """
    def __init__(self, ensemble ):
        self.inputformat = 'Towhee'
        self.randomseed = int(random.random()*10000000)
        self.irandom_luxlevel = 3
        self.random_allow_restart = 'T'
        self.ensemble = ensemble
        self.temperature = 300.0
        self.pressure = 95.0
        self.nmolty = 0
        self.nmolectyp = []
        self.chempot = 0.0
        self.numboxes = 1
        self.stepstyle = 'moves'
        self.nstep = 0
        self.printfreq = 0
        self.blocksize = 0
        self.moviefreq = 0
        self.backupfreq = 0
        self.runoutput = 'full'
        self.pdb_output_freq = 0
        self.louthist = 'T'
        self.hist_lable = 1
        self.hist_suffix = 'd'
        self.his_nequil = 10
        self.histcalcfreq = 10
        self.pressurefreq = 10
        self.trmaxdispfreq = 1000
        self.volmaxdispfreq = 1000
        self.potentialstytle = 'internal'
        self.ffnumber = 1
        self.filename = ''
        self.classical_potential = 'Lennard-Jones'
        self.calssica_mixrule = 'Lorentz-Berthelot'
        self.lshift = '.false.'
        self.ltailc = '.true.'
        self.rmin = 1.0
        self.rcut = 10.0
        self.rcutin = 10.0
        self.electrostatic_form = 'coulomb'
        self.coulombstyle = 'ewald_fixed_kmax'
        self.kalp = 5.6
        self.kmax = 5
        self.dielect = 1.0
        self.linit = 'T'
        self.initboxtype = ['dimensions',]
        self.initstyle = ['full cbmc', ]
        self.initlattice = ['simple cubic']
        self.initmol = 1
        self.inin = [0, 0, 0]
        self.hmatrix = [[0,0,0],]
        # mc moves
        self.pmuvtcbswap = 0.0
        self.pmuvtcbmt = [1.0]
        self.pmcb = 0.0
        self.pmcbmt = [1.0]
        self.pmall = [0.0]
        pmtracm = 0
        self.pmtcmt = [1.0]
        self.rmtrac = [0.5]
        self.tatrac = [0.5]
        self.pmrotate = 0.0
        self.pmromt = [1.0]
        self.rmrot = [0.05]
        self.tarot = [0.5]
        self.icbmc_formulation = 'Martin and Frischknecht 2006'
        self.cbmc_setting_style = 'Martin and Frischknecht'
        self.input_style = 'input_style'
        self.nunit = 0
        self.nmaxcbmc = 0
        self.lpdbnames = 'F'
        self.forcefield = 'TEAM'
        self.charge_assignment = 'manual'
        # system
        self.sys = System()

    def writeInput(self,):
        # check if there are molecules
        if len(self.sys.mols) == 0:
            print 'Warning: No molecules are in the system!'
            return 0
        if len(self.sys.atoms) == 0:
            print 'Warning: No atoms are in the system!'
            return 0



    

a = Towhee('uvt')
a.writeInput()

