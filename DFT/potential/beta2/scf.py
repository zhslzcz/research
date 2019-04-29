#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 10:26:33 2018

@author: clian
"""

from ase.io import read, write
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.emt import EMT
import numpy as np

from gpaw import GPAW
from math import sqrt
import numpy as np
from ase import Atoms
from ase.parallel import paropen
from gpaw import GPAW, FermiDirac, PW
from ase.io import write, read

from ase.io.trajectory import Trajectory
atoms = read('beta2.relaxed.vasp')
calc = GPAW(mode='lcao', basis='dzp',#mode=PW(350),
            kpts=(2, 1, 1), xc='PBE', occupations=FermiDirac(0.05))
atoms.set_calculator(calc)
atoms.get_potential_energy()
calc.write('scf.gpw', 'all')
