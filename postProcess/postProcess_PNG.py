#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 16:26:54 2021

@author: yunjaeh
"""
import numpy as np
import matplotlib.pyplot as plt

class fCase():
    def __init__(self):
        self.X=[]
        self.Y=[]
        self.Z=[]
        self.U=[]
        self.V=[]
        self.W=[]
        self.P=[]
        
    def readData(self,fPath):
        data_temp = np.loadtxt(fPath+'U_avg.raw_values.dat')
        self.X = data_temp[:,4]
        self.Y = data_temp[:,5]
        self.Z = data_temp[:,6]
        self.U=data_temp[:,-1]
        data_temp = np.loadtxt(fPath+'V_avg.raw_values.dat')
        self.V=data_temp[:,-1]
        data_temp = np.loadtxt(fPath+'W_avg.raw_values.dat')
        self.W=data_temp[:,-1]
        data_temp = np.loadtxt(fPath+'P_avg.raw_values.dat')
        self.P=data_temp[:,-1]
        
        

#%%
fPath='../LES_ventilation/E1_5/post_avg/'
E1_5 = fCase()
E1_5.readData(fPath)

fPath='../LES_ventilation/E1_5.base_inflow/post_avg/'
E1_5_inflow = fCase()
E1_5_inflow.readData(fPath)


#%%        
from scipy.interpolate import griddata

unique_X = np.unique(E1_5.X)
unique_Y = np.unique(E1_5.Y)
# len_X = len(unique_X)
# len_Y = len(unique_Y)

XY=np.transpose([E1_5.X, E1_5.Y])

delta_x=0.005
grid_x, grid_y = np.mgrid[0.4:0.7:delta_x, 0:0.2:delta_x ]
grid_u=griddata(XY,E1_5.U,(grid_x,grid_y))
grid_v=griddata(XY,E1_5.V,(grid_x,grid_y))

grid_u2=griddata(XY,E1_5_inflow.U,(grid_x,grid_y))
grid_v2=griddata(XY,E1_5_inflow.V,(grid_x,grid_y))


#%%
D=0.1
U_ref=6.6
E1_5.mag = np.sqrt(E1_5.U**2.0  + E1_5.V**2.0)
E1_5_inflow.mag = np.sqrt(E1_5_inflow.U**2.0  + E1_5_inflow.V**2.0)

fig, axes = plt.subplots(nrows=1, ncols=2,figsize=(12,4))

im0 =axes[0].scatter((E1_5.X-0.5)/D,E1_5.Y/D,10, E1_5_inflow.mag/U_ref)
axes[0].quiver((grid_x-0.5)/D,grid_y/D,grid_u2,grid_v2, scale=100)

im1 =axes[1].scatter((E1_5.X-0.5)/D, E1_5.Y/D,10,E1_5.mag/U_ref)
axes[1].quiver((grid_x-0.5)/D, grid_y/D,grid_u,grid_v, scale=100)

im0.set_cmap(cmap=plt.get_cmap('jet'))
im1.set_cmap(cmap=plt.get_cmap('jet'))
im0.set_clim(0,1.2)
im1.set_clim(0,1.2)

axes[0].set_title('LES, base inflow')
axes[1].set_title('LES, optimized inflow')
for i in range(0,2):
    axes[i].set(xlim=(-0.05/D,0.15/D), xticks=[-0.5,0,0.5,1,1.5], \
                ylim=(0,0.15/D), yticks=[0,0.5,1,1.5], \
                xlabel='x/D', ylabel='y/D')

fig.colorbar(im0, ax=axes, orientation='vertical')
fig.savefig('E1_5_inflow_velocity_magnitude.png')

#%%
fig, axes = plt.subplots(nrows=1, ncols=2,figsize=(12,4))
axes[0].quiver((grid_x-0.5)/D,grid_y/D,grid_u2,grid_v2, scale=100)



#%%
fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(6,4))
# axes.quiver((grid_x-0.5)/D,grid_y/D,grid_u2,grid_v2, color='k', scale=200)
axes.quiver((grid_x-0.5)/D, grid_y/D,grid_u,grid_v,color='r',scale=200)
axes.quiver((grid_x-0.5)/D,grid_y/D,grid_u2,grid_v2, color='k', scale=200)
axes.set(xlim=(-0.05/D,0.15/D), xticks=[-0.5,0,0.5,1,1.5], \
                ylim=(0,0.15/D), yticks=[0,0.5,1,1.5], \
                xlabel='x/D', ylabel='y/D')


#%%
data_ref = np.loadtxt('../ReferenceData/E1.csv',skiprows=2,delimiter=',')
idx_ref = data_ref[:,0] < 0.0
data_ref[idx_ref,0] = data_ref[idx_ref,0]*0.25/0.35


#%%
Uref=6.6
D=0.1
# idx = (E1_5.Y == 0.0396875) 
idx1 = (E1_5.Y == 0.0403125) 
idx2 = (E1_5_inflow.Y == 0.0403125) 

plt.plot(data_ref[:,0],data_ref[:,1],'k.',markersize=10)
plt.plot((E1_5_inflow.X[idx1]-0.5)/D, E1_5_inflow.U[idx1]/Uref,'r',linewidth=2)
plt.plot((E1_5.X[idx1]-0.5)/D, E1_5.U[idx1]/Uref,'b',linewidth=2)

plt.legend(['PIV','LES: base inflow','LES: optimized inflow'])
plt.xlim([-0.25, 1.25])
plt.ylim([0, 1])
plt.xticks([-0.25,0,1,1.25])
plt.xlabel('x/D')
plt.ylabel('U/U$_{ref}$')
plt.grid()

# plt.savefig('../Results/E1_5_validation.png')




