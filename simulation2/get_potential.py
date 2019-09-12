from ase.io.cube import read_cube_data
import pandas as pd
import numpy as np
# from ase import Atoms
# from ase.visualize import view
import matplotlib.pyplot as plt

#constants:

dx,dy,dz = 8.54/44,17.08/84,25/124
cellX,cellY,cellZ=8.54, 17.08, 25.


class DFTpotential(object):
	"""process the dft file"""
	def __init__(self, data_path,fermishift): #'../zeta_A_potential.cube' or '../beta_2_potential.cube'
		self.path = data_path
		self.potential,self.atoms=read_cube_data(self.path)
		self.offset = self.potential.mean(axis=0).mean(axis=0).max()+ fermishift#potential at infinite far
		self.atomxyz = self.atoms.get_positions()
		self.symbols = self.atoms.get_chemical_symbols()
		self.atomdf = pd.concat([pd.DataFrame(self.atomxyz,columns=['x', 'y', 'z']),pd.DataFrame(self.symbols,columns=['atom'])], axis=1)

	def normXY(self,x,y):
		#convert XYZ to within boundary first
	    while  x<0: 
	        x+=cellX 
	    while  x>cellX: 
	        x-=cellX
	    while  y<0: 
	        y+=cellY 
	    while  y>cellY: 
	        y-=cellY

	    return x,y

	def potentialAtXYZ(self,x,y,z):
		'''determine the potential at any xyz'''
		x,y = self.normXY(x,y)
		X,Y,Z = int(x/dx),int(y/dy),int(z/dz)
		#print(X,Y,Z)
		return round(self.potential[X][Y][Z]-self.offset,3) 

	def potential_vs_z(self,x,y,z):
		'''get potential above (x,y,z)'''
		x,y = self.normXY(x,y)
		X,Y,Z = int(x/dx),int(y/dy),int(z/dz)
		return np.round(self.potential[X][Y][Z:]-self.offset,3) #the surface has largest z

	def avg_potential_vs_z(self,x,y,z,r):
		'''get average potential in radius r'''
		R=int(r/dx) # dx ,dy are both around 0.2
		lines=[]
		for delta_x in range(-R,R+1):
			for delta_y in range(-R,R+1):
				nx,ny = x+delta_x*dx,y+delta_y*dy
				lines.append(self.potential_vs_z(nx,ny,z))
		return np.round(np.mean(lines,axis=0),3)

	def plot_potential(self,potential):

		plt.plot(potential)




# zeta = DFTpotential('../zeta_A_potential.cube')
# print(zeta.potential_vs_z(5,5,10))

# p = zeta.avg_potential_vs_z(5,5,10,1)
# print(p.shape)
# plt.figure()
# zeta.plot_potential(p)
# plt.show()




