import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class LEISsim(object):
    """calculates the scattering probability from LEIS simulation result"""
    def __init__(self, LEIS_file_path, dft_structure_path): #'../dftzetaA0.txt','../zetaA0.txt'
        self.LEIS = LEIS_file_path
        self.dft = dft_structure_path
        self.atoms = pd.read_table(self.dft,
                                    delim_whitespace=True,
                                    names = ('x','y','z','Z','m','f')
                                    )
        self.LEIStable = pd.read_table(self.LEIS)
        self.toplist = self.unblockedatoms()
        self.X = np.array(self.LEIStable['rx'])
        self.Y = np.array(self.LEIStable['ry'])
        self.bx = max(self.X)-min(self.X)
        self.by = max(self.Y)-min(self.Y)
        self.scatterWeight()


    def scatterWeight(self):
        counts = np.zeros(self.atoms.shape[0])
    
        topx=np.array([self.atoms.loc[a]['x'] for a in self.toplist])
        topy=np.array([self.atoms.loc[a]['y'] for a in self.toplist])

        for x,y in zip(self.X,self.Y):
            lx = np.abs(topx-x)
            ly = np.abs(topy-y)
            ls = np.minimum(lx,self.bx-lx)**2+np.minimum(ly,self.by-ly)**2
            index_min = self.toplist[np.argmin(ls)]
            counts[index_min] += 1
        self.atoms['count'] = counts
        self.atoms['sum'] = self.atoms.groupby('Z')['count'].transform('sum')
        self.atoms['perc'] = (self.atoms['count']/self.atoms['sum']).round(3)

    def unblockedatoms(self,block_radius=0.4):#find the top unblocked atoms
        toplist = []    
        for index, row in self.atoms.iterrows():
            k=0
            for a in toplist:
                rowa = self.atoms.loc[a]
                if (rowa['x']-row['x'])**2+(rowa['y']-row['y'])**2<block_radius**2:
                    if rowa['z']>row['z']:
                        k=1
                        break
                    else:
                        toplist.remove(a)
            if k==0:
                toplist.append(index)
        return toplist

    def AtomXYZ(self,boundary=True,limit=0.02):
        topIn = self.atoms.loc[self.toplist].query('Z ==49')
        topAs = self.atoms.loc[self.toplist].query('Z ==33')
        In_X = list(topIn['x'])
        In_Y = list(topIn['y'])
        As_X = list(topAs['x'])
        As_Y = list(topAs['y'])
        
        size_In = [240 if z>-2 else 80 for z in topIn['z']]
        size_As = [240 if z>-2 else 80 for z in topAs['z']]
        
        
        minx=min(self.X)
        maxx=max(self.X)
        miny=min(self.Y)
        maxy=max(self.Y)

        if boundary:
        
            tempIn = zip(list(In_X),list(In_Y))
            limit=0.02
            for i,(x,y) in enumerate(tempIn):
                if x-minx<limit:
                    In_X.append(x+self.bx)
                    In_Y.append(y)
                    size_In.append(size_In[i])
                    
                    if y-miny<limit:
                        In_X.append(x+self.bx)
                        In_Y.append(y+self.by)
                        size_In.append(size_In[i])
                
                if maxx-x<limit:
                    In_X.append(x-self.bx)
                    In_Y.append(y)
                    size_In.append(size_In[i])
                
                if y-miny<limit:
                    In_X.append(x)
                    In_Y.append(y+self.by)
                    size_In.append(size_In[i])
                
                if maxy-y<limit:
                    In_X.append(x)
                    In_Y.append(y-self.by)
                    size_In.append(size_In[i]) 
            
            tempAs = zip(list(As_X),list(As_Y))
            for i,(x,y) in enumerate(tempAs):
                if x-minx<limit:
                    As_X.append(x+self.bx)
                    As_Y.append(y)
                    size_As.append(size_As[i])
                    
                    if y-miny<limit:
                        As_X.append(x+self.bx)
                        As_Y.append(y+self.by)
                        size_As.append(size_As[i])
                    
                    if maxy-y<limit:
                        As_X.append(x+self.bx)
                        As_Y.append(y-self.by)
                        size_As.append(size_As[i]) 

                if maxx-x<limit:
                    As_X.append(x-self.bx)
                    As_Y.append(y)
                    size_As.append(size_As[i])
                    
                
                if y-miny<limit:
                    As_X.append(x)
                    As_Y.append(y+self.by)
                    size_As.append(size_As[i])
                
                if maxy-y<limit:
                    As_X.append(x)
                    As_Y.append(y-self.by)
                    size_As.append(size_As[i])


        #print(XYZ)
        return In_X,In_Y,size_In,As_X,As_Y,size_As

    def plottopAtoms(self,**kwarg):
        
        plt.figure(figsize=(self.bx/2,self.by/2))
        plt.xticks([])
        plt.yticks([])

        plt.scatter(self.X,self.Y,marker='.',alpha=0.5,label='scatter sites')
        
        In_X,In_Y,size_In,As_X,As_Y,size_As = self.AtomXYZ()

        plt.scatter(In_X,In_Y,marker='o',s=size_In,color='orange',label='In')
        plt.scatter(As_X,As_Y,marker='o',s=size_As,color='g',label='As')
        #    plt.axvline(0.00000,linestyle='--',color='grey')
        plt.axhline(0.00000,linestyle='--',color='grey')
        plt.legend(bbox_to_anchor=(1.05, 1),fontsize=15)
        plt.show()

    def write_to_csv(self,name):
        self.atoms.to_csv(name)

for [LEIS_file_path,dft_structure_path,name] in (['../dftzetaA0.txt','../zetaA0.txt','zetaA0'],
                                                 ['../dftzetaA90.txt','../zetaA90.txt','zetaA90'],
                                                 ['../dftbeta0.txt','../beta0.txt','beta0'],
                                                 ['../dftbeta90.txt','../beta90.txt','beta90']):

    zetaA0 = LEISsim(LEIS_file_path,dft_structure_path)
    zetaA0.write_to_csv('./weights/'+name)
    #print(zetaA0.atoms)





    
            
    
    
    
