import numpy as np
from scipy.constants import value, Planck
import ase.io
from panther.io import read_vasp_hessian
from panther.vibrations import harmonic_vibrational_analysis
from panther.thermochemistry import Thermochemistry

meoh = ase.io.read('POSCAR')

hessian = read_vasp_hessian('OUTCAR', symmetrize=True, convert_to_au=True)
np.save('hessian', hessian)

frequencies, normal_modes = harmonic_vibrational_analysis(hessian, meoh, proj_translations=True, proj_rotations=True, ascomplex=False)

vibenergies = Planck * frequencies.real * value('hartree-hertz relationship')
print(vibenergies)
vibenergies = vibenergies[vibenergies > 0.0]

thermo = Thermochemistry(vibenergies, meoh, phase='gas', pointgroup='Cs')
thermo.summary(T=273.15, p=0.1)
