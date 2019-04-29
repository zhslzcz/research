#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 21:29:16 2018

@author: clian
"""

import numpy as np
from matplotlib import pyplot as plt
from ase.io import write
from ase.dft.stm import STM
from gpaw import GPAW


calc = GPAW('scf.gpw', txt=None)
atoms = calc.get_atoms()
stm = STM(atoms)
bias = -0.5
stm.calculate_ldos(bias=bias)
write('ldos_bias='+ str(bias) +'.cube',atoms, data=stm.ldos)

# 2d plot
# c = stm.get_averaged_current(bias=1.0, z=19)
# x, y, h = stm.scan(bias, c, repeat=(3, 2))

# import matplotlib.pyplot as plt

# plt.gca(aspect='equal')
# plt.contourf(x, y, h, 40)
# plt.colorbar()
# plt.show()