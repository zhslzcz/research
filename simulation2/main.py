import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import get_potential as gp 
import freeze_model as fm

# tunable parameters:

fzd = {'In': fm.freezingDistance('In'),'As': fm.freezingDistance('As')} #{2.83,2.94}
fermishift = {'zetaA': 3, 'beta':3.8}
R_image_plane = 5
r_potential_avg = 3
fz_range = np.linspace(1,4,20)
image_plane_offset = 0


class Table_Update(object):
	"""update panda table """
	def __init__(self,structure,angle,R=R_image_plane,r=r_potential_avg,fz_range=fz_range): 
		if structure in ['zetaA','beta']:
			self.structure = structure
		else:
			raise 'wrong structure name'
		if angle in [0,90]:
			self.angle=angle
		else:
			raise 'wrong angle'
		self.r = r # the radius of area that averages the potential around
		self.R = R # the image plane search radius for contour model
		self.fz_range=fz_range
		self.weight = pd.read_csv('./weights/'+structure+str(angle))
		self.potential = self.getpotential(fermishift[self.structure])

		self.table = pd.concat([self.potential.atomdf,self.weight[['count','perc']]],axis=1)
		self.getpotentialatfz()
		self.getimageZ()
		self.report = self.getreport()

	def getimageZ(self):
		imageZs = []
		for x,y in self.table[['x','y']].values:
			imageZ = fm.outmostZ(self.potential.atomxyz,x,y,self.R)+1
			imageZs.append(imageZ)

		self.table['image_Z'] = imageZs

	def getpotential(self,shift):
		if self.structure =='zetaA':
			data_path = '../zeta_A_potential.cube'
		else:
			data_path = '../beta_2_potential.cube'

		return gp.DFTpotential(data_path,shift)

	def getpotenialabove(self):
		p_vs_z=[]
		for x,y,z in self.table[['x','y','z']].values:
			p_vs_z.append(self.potential.avg_potential_vs_z(x,y,z,self.r))
		self.table['p_vs_z']= p_vs_z

	def getNFabove(self):
		NFs=[]
		for x,y,z in self.table[['x','y','z']].values:
			NF = fm.contourModel(self.potential,self.r,x,y,z, squence=True)
			NFs.append(NF)
		self.table['NF_above'] = NFs

	def NFatXYZ(self,x,y,z):
		NF = fm.contourModel(self.potential,self.R,x,y,z,self.r,image_plane_offset)
		return NF

	def getNFatFreezingDistance(self):
		NFs = []
		for x,y,z,spicie,image_Z in self.table[['x','y','z','atom','image_Z']].values:
			fz = fzd[spicie]
			NFs.append(self.NFatXYZ(x,y,image_Z+fz)) #flatmodel
			#NFs.append(self.NFatXYZ(x,y,z+fz)) # contourmodel

		self.table['NF_frzpoint'] = NFs

	def getpotentialatfz(self):
		ps = []
		for x,y,z,spicie in self.table[['x','y','z','atom']].values:
			fz = fzd[spicie]
			ps.append(self.potential.potentialAtXYZ(x,y,z+fz))
		self.table['p_frzpoint'] = ps

	def getweightedNF(self):
		if 'NF_frzpoint' not in self.table.columns:
			self.getNFatFreezingDistance()
		self.table['weight_NF_frzpoint'] = self.table['NF_frzpoint']*self.table['perc']

	def get_NF_vs_fz(self):
		'''calculate NF at differnt freeze distance'''
		NF_vs_fz = []
		for x,y,z,spicie,image_Z in self.table[['x','y','z','atom','image_Z']].values:
			NFs = []
			for fz in self.fz_range:
				#NFs.append(self.NFatXYZ(x,y,z+fz)) #contour model
				NFs.append(self.NFatXYZ(x,y,image_Z+fz)) # flatmodel
			NF_vs_fz.append(NFs)
		self.table['NF_vs_fz'] = NF_vs_fz

	def getweightedNF_vs_fz(self):
		if 'NF_vs_fz' not in self.table.columns:
			self.get_NF_vs_fz()
		weighted=[]
		for perc,value in self.table[['perc','NF_vs_fz']].values:
			weighted.append(perc*np.array(value))
		self.table['weight_NF_vs_fz'] = weighted

	def getreport(self):
		self.getNFatFreezingDistance()
		self.getweightedNF()
		self.getweightedNF_vs_fz()
		grouped = self.table.groupby(['atom'])
		report = grouped[['perc','weight_NF_frzpoint']].agg('sum')
		temp = grouped[['weight_NF_vs_fz']].agg(lambda x: list(x))
		weighted=[]
		for value in temp['weight_NF_vs_fz'].values:
			weighted.append(np.sum(value,axis=0).round(1))
		report['weight_NF_vs_fz'] = weighted

		return report

	def plot_NF_vs_fz(self):
		linestyle = {'In':'--','As':'-'}
		for atom,row in self.report.iterrows():
			NFs = row['weight_NF_vs_fz']
			plt.plot(self.fz_range,NFs,linestyle[atom],label=atom+str(self.angle))
		

expvalues={'zetaA':{0:{'As':57,'In':48,'As_e':6,'In_e':4},90:{'As':65,'In':40,'As_e':12,'In_e':3}},
	   'beta':{0:{'As':66,'In':34,'As_e':14,'In_e':3},90:{'As':50,'In':41,'As_e':6,'In_e':5}}}

def exp_plot(data):
	structure,angle = data.structure,data.angle
	expNF= expvalues[structure]
	for a,c,l in [('As','b','-'),('In','g','--')]:
		plt.axhline(y=expNF[angle][a], color=c, linestyle=l)
		plt.axhspan(expNF[angle][a]-expNF[angle][a+'_e'],expNF[angle][a]+expNF[angle][a+'_e'],color=c, linestyle=l, alpha = 0.1)
		plt.axvline(x=fzd[a],color='k', linestyle=l)

for s in ['zetaA','beta']:
	for a in [0,90]:
		print(s,a)
		m=Table_Update(s,a)
		m.getpotentialatfz()
		m.getNFatFreezingDistance()
		print(m.getreport())
		m.table.drop(['NF_vs_fz','weight_NF_vs_fz'],axis=1,inplace=True)
		# print(m.report['weight_NF_frzpoint'].values)
		# print(m.table)
		fig,ax = plt.subplots()
		ax.figsize=(10,10)
		plt.axis([m.fz_range[0], m.fz_range[-1], 0, 100])
		m.plot_NF_vs_fz()
		exp_plot(m)
		plt.legend(loc='best')
		ax.set_title(s+str(a))
		
plt.show()


# exp = [[57,48,6,4],[65,40,12,3],[66,34,14,3],[50,41,6,5]]




# fz_range = np.linspace(1,4,20)
# for zetashift in np.linspace(-2,1,5):
# 	for betashift in np.linspace(-2,1,5):
# 		print(zetashift,betashift)
# 		fermishift = {'zetaA': zetashift, 'beta':betashift}
# 		tables=[]
# 		for s in ['zetaA','beta']:
# 			for a in [0,90]:
# 				tables.append(Table_Update(s,a))

# 		for k in np.linspace(0.9,1.2,5):
# 			fzd = {'In': 2.83*k,'As': 2.94*k}	
# 			for R_image_plane in np.linspace(0,10,3):
# 				for r_potential_avg in np.linspace(0,10,5):
# 					for image_plane_offset in np.linspace(0,2,5):
# 						match=0
# 						for i,table in enumerate(tables):
# 							report = table.getreport()
# 							expAs,expIn,expAs_e,expIn_e=exp[i]
# 							As,In = report['weight_NF_frzpoint'].values
# 							if expAs-expAs_e<As<expAs+expAs_e and expIn-expIn_e<In<expIn+expIn_e:
# 								match+=1
# 						if match>=2:
# 							print('match found',match,k,fermishift,R_image_plane,r_potential_avg,image_plane_offset)







