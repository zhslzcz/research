from __future__ import print_function
import numpy as np
import pylab as plt

from gpaw import GPAW
from ase.io import write

# Read in the 5-layer slab:
calc = GPAW('scf.gpw', txt=None)
slab = calc.get_atoms()

# Get the height of the unit cell:
L = slab.get_cell()[2, 2]

# Get the effective potential on a 3D grid:
v = calc.get_electrostatic_potential()


write('potential.cube', slab, data=v)