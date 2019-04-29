#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 21:29:16 2018

@author: clian
"""

import numpy as np
from matplotlib import pyplot as plt
from py3ramids.plot.PlotUtility import scanFolder
import py3ramids.plot.setting as ma
from ase.io import write
def getData(index, folder):
    from gpaw import GPAW
    import os
    if not os.path.exists('1DPot.npz'):   
        # Read in the 5-layer slab:
        calc = GPAW('scf.gpw', txt=None)
        slab = calc.get_atoms()
        efermi = calc.get_fermi_level()
        # Get the height of the unit cell:
        
        # Get the effective potential on a 3D grid:
        v = calc.get_effective_potential()
#        for 
#        np.save('1DPot', slab=slab, v=v, efermi=efermi)
#    else:
#        slab, v, efermi = np.load('1DPot.npz')
#    
    return slab, v, efermi
    
# start the plotting code 
def plot(iv, ax):
    data = scanFolder(getData)
    labels=[r'$\beta_2$',r'$\zeta_A$']
    for iline, (slab, v, efermi) in enumerate(data):
        L = slab.get_cell()[2, 2]
        zatoms = slab.positions[:,2]
        print(zatoms.max()-zatoms.min())
        nx, ny, nz = v.shape
        print(v.shape)
        vz = v.mean(axis=0).mean(axis=0)
        z = np.linspace(0, L, nz, endpoint=False)
        write('%s_potential.cube'%labels[iline],slab, data=v)
        ax.plot(z-zatoms.min(), vz-vz[0], label=labels[iline])
        
    label = ''
    zSurf = 0
    kargs=ma.getPropertyFromPosition(xlabel='',
                                     ylabel=r'',
                                     vline=[zSurf,zSurf+3],
                                     title='')
      
      
    ma.setProperty(ax,**kargs)
    #ma.add_label(iv, ax)
    return kargs
# end the plotting code 


if __name__ == '__main__':
    fig, axs = plt.subplots(1,1,sharex=False,sharey=False,figsize=(8,6))
    plot(0,axs)
    plt.tight_layout()
    fig.align_ylabels(axs)
    SaveName = __file__.split('/')[-1].split('.')[0]
    if True:
      for save_type in ['.pdf','.png']:
        filename = SaveName + save_type
        plt.savefig(filename,dpi=600)



        