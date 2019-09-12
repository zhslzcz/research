import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import get_potential as gp 




#constants:
M = {'Na': 22.99,'In':114.82,'As':74.92}
vmax= 2.6

def freezingDistance(spicie,E=3000,theta=0,sc=125):
	m1 = M['Na']
	m2 = M[spicie]
	E = E*((m1*np.cos(sc/180*np.pi)+np.sqrt(m2**2-(m1*np.sin(sc/180*np.pi))**2))/(m1+m2))**2
	#print (E)
	v = np.sqrt(8389769*E)*np.cos(theta/180*np.pi)
	#print(v/2187691)
	d = 1/0.86*np.log(2*2.23/(0.86*v/2187691))*0.53
	return round(d,2)

def energyShift(d_to_image):
	'''calculate energy shifting using distance to image plane '''
	if d_to_image<0:
		return 0
	else:
		return 27.2/np.sqrt((27.2*27.2)/(vmax*vmax)+16*(d_to_image/0.53)**2) 

def broadening(d_to_image):
	'''calculating broadening of the s level'''
	if d_to_image<0:
		return 27.2*0.04
	else:
		return 27.2*2.23/(np.exp(4*0.86*d_to_image/0.53)+(2.23/0.04)**4-1)**0.25

def potentialToNF(potential,d_to_image):

	Eshift = energyShift(d_to_image)
	Ebroad = broadening(d_to_image)
	Es = -5.14+Eshift
	# guassion integral
	sigma = Ebroad/1.38629
	NF = 0.5+0.5*np.erf((potential-Es)/(1.414*sigma))
	return np.round(NF*100,1) #NF in %

def outmostZ(atomxyz,X,Y,R):
	'''calculate the outmost atom_z at (X,Y) in radius r'''
	Zmax=0
	for x,y,z in atomxyz:
		if z>Zmax and min(abs(x-X),gp.cellX-abs(x-X))**2+min(abs(y-Y),gp.cellY-abs(y-Y))**2<=R**2:
			Zmax = z
	return Zmax

def contourModel(DFT,R,x,y,z,r,image_plane_offset=1,squence=False):
	'''use RCT model to determine NF at (x,y,z)'''
	image_z = outmostZ(DFT.atomxyz,x,y,R)+image_plane_offset
	d_to_image = z-image_z	
	if squence:
		potential = DFT.avg_potential_vs_z(x,y,z,r)
	else:
		potential = DFT.potentialAtXYZ(x,y,z)
	NF = potentialToNF(potential,d_to_image)
	return NF

	
# zeta = gp.DFTpotential('../zeta_A_potential.cube')
# plt.figure()
# for x,y,z in zeta.atomxyz:
# 	NF = contourModel(zeta,2,x,y,z, squence=True)
# 	z = z+np.arange(len(NF))*gp.dz
# 	plt.plot(z,NF)
# plt.show()





	
		