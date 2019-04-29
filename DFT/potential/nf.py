#!/usr/bin/env python
# coding: utf-8


from ase.io.cube import read_cube_data
import pandas as pd
import numpy as np
from ase import Atoms
from ase.visualize import view
import math
import matplotlib.pyplot as plt
                      

#calculate freezing distancce of 3000eV Na scattered from InAs at 125º
#the projectile comes out in surface normal direction
def freezingDistance(E):
    #NA=6.02*(10^23)   #Avogadora 
    #m=22.989769/NA/1000 #Na mass in kg
    #E=E*1.60217662*(10^-19) #Energy in Joul
    #v=math.sqrt(2*E/m)=math.sqrt(2*E*1.602/22.99*6.02*10000000)
    v=math.sqrt(8389769*E)
    d=1/0.86*math.log(2*2.23/(0.86*v/2187691))*0.53
    return round(d,2)


# determine the potential at freezing distance
def potential(row):
    X,Y = int(round(row['x']/dx)),int(round(row['y']/dy))
    zf = freeze_In if row['atom'] == 'In' else freeze_As         
    Z = int(round((row['z'] + zf) / dz))
    
    return data[X][Y][Z]- offset  #not sure if this is right


# determine the potential above each atom
def potential_vs_z(row,Zrange):
    
    X,Y,Z = row['x'],row['y'],row['z']
    p_vs_z=[]
    for delta_z in Zrange: 
        p_vs_z+=[round(data[int(X/dx)][int(Y/dy)][int((Z+delta_z)/dz)]- offset,2)]
    return p_vs_z 


# determine the average potential above each atom
def avg_potential_vs_z(row,Zrange):
    X,Y,Z = row['x'],row['y'],row['z']
    avg_p_z=[]
    r=1    #radius of the avarage area
    XYpairs=[]
    for delta_x in np.linspace(-r,r,9):
        for delta_y in np.linspace(-r, r, 9):
            if np.abs(delta_x+delta_y*1.j)<=r:

                PX,PY=X+delta_x,Y+delta_y

                #deal with out of boundary condition
                while  PX<0: 
                    PX+=cellX 
                while  PX>cellX: 
                    PX-=cellX
                while  PY<0: 
                    PY+=cellY 
                while  PY>cellY: 
                    PY-=cellY
                
                XYpairs+=[(PX,PY)]
                
    for delta_z in Zrange: 
        p_tot = []
        for PX,PY in XYpairs:
            p_tot.append(data[int(PX/dx)][int(PY/dy)][int((Z+delta_z)/dz)])
        
        avg_p_z+=[round(np.mean(p_tot) - offset,2)]

    return avg_p_z 


#use RCT model to determine NF scattered from each site
def NF(row):
    zf = freeze_In if row['atom'] == 'In' else freeze_As  
    
    #calculate s level after shifting 
    dE = 27.2/math.sqrt((27.2*27.2)/(2.6*2.6) + 16*(zf/0.53)**2)
    Es = -5.14+dE
    
    #calculating broadening of the s level
    broad = 27.2*2.23/math.sqrt(math.sqrt(math.exp(4*0.86*zf/0.53)+(2.23/0.04)**4-1))
    
    #guassion integral
    sigma = broad/1.38629
    NF = 0.5+0.5*math.erf((row['dz_potential']-Es)/(1.414*sigma))
    
    return round(NF*100,1) #NF in %


# some constants

freeze_As=freezingDistance(1100) #freezing distance from As
freeze_In=freezingDistance(1582) #freezing distance from In
#freeze_As=1.6 #freezing distance from As
#freeze_In=1.7 #freezing distance from In
print(freeze_As,freeze_In)

dx,dy,dz = 8.54/44,17.08/84,25/124
cellX,cellY,cellZ=8.54, 17.08, 25.


#dataname='zeta_A_potential.cube'
dataname='beta_2_potential.cube'

data,atoms=read_cube_data(dataname)

if dataname=='zeta_A_potential.cube':
    visible_90_As = [0,1,2,3,12,13,14,15]
    visible_90_In = [4,5,6,7,8,9,10,11,19,23]
    visible_0_As = [0,1,2,3,12,13,14,15]
    visible_0_In = [4,5,6,7,8,9,10,11]
    efermi = - 0.2086
    wf = 5.24
else:
    visible_90_As = [6,7,8,9,15,17,44,45,46,47]
    visible_90_In = [10,11]
    visible_0_As = [6,7,8,9,46,47]
    visible_0_In = [10,11,12,13,42,43,48,49,50,51]
    efermi = - 0.322
    wf = 5.38

offset = data.mean(axis=0).mean(axis=0).max()  #potential at infinite far
    
atomxyz=atoms.get_positions()
symbols=atoms.get_chemical_symbols()

atomdf = pd.concat([pd.DataFrame(atomxyz,columns=['x', 'y', 'z']),
                    pd.DataFrame(symbols,columns=['atom'])],
                   axis=1)
# print(atoms)
# print(data.shape)
# view(atoms)
# plt.imshow(data[20])
# plt.show()


atomdf['dz_potential']=atomdf.apply(potential,axis=1)
atomdf['NF']=atomdf.apply(NF,axis=1)


atomdf.loc[visible_90_As]

atomdf.loc[visible_90_In]


Zrange=np.linspace(0,4,20)
atomdf['p_vs_z']=atomdf.apply(potential_vs_z,Zrange=Zrange,axis=1)



for i in visible_90_As:
    plt.plot(Zrange,atomdf.loc[i]['p_vs_z'],'-',label='As'+str(i))
for i in visible_90_In:
    plt.plot(Zrange,atomdf.loc[i]['p_vs_z'],'*', label='In'+str(i))
# for i in visible_90_As:
#     plt.plot(Zrange,atomdf.loc[i]['p_vs_z'],'--',label='As'+str(i))
# for i in visible_90_In:
#     plt.plot(Zrange,atomdf.loc[i]['p_vs_z'],'+', label='In'+str(i))
#plt.legend(loc='best')
#plt.savefig(dataname[:4]+'_p_vs_z')


atomdf['avg_p_z']=atomdf.apply(avg_potential_vs_z,Zrange=Zrange,axis=1)


for i in visible_90_As:
    plt.plot(Zrange,atomdf.loc[i]['avg_p_z'],'-',label='As'+str(i))
for i in visible_90_In:
    plt.plot(Zrange,atomdf.loc[i]['avg_p_z'],'*', label='In'+str(i))
#plt.legend(loc='best')
#plt.savefig(dataname[:4]+'avg_p_vs_z')




