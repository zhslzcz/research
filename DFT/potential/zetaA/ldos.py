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
stm = STM(slab)
bias = 0.5
stm.calculate_ldos(bias=bias)
write('ldos_bias='+ bias +'.cube',atoms, data=stm.ldos)