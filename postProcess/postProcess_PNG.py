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

fPath='../LES_ventilation/E1_5.fine/post_avg/'
E1_5_fine = fCase()
E1_5_fine.readData(fPath)

fPath='../LES_ventilation/E1_5.coarse/post_avg/'
E1_5_coarse = fCase()
E1_5_coarse.readData(fPath)


#%%        
from scipy.interpolate import griddata

unique_X = np.unique(E1_5.X)
unique_Y = np.unique(E1_5.Y)

XY=np.transpose([E1_5.X, E1_5.Y])

delta_x=0.005
grid_x, grid_y = np.mgrid[0.4:0.7:delta_x, 0.001:0.2:delta_x ]
grid_u=griddata(XY,E1_5.U,(grid_x,grid_y))
grid_v=griddata(XY,E1_5.V,(grid_x,grid_y))
grid_w=griddata(XY,E1_5.W,(grid_x,grid_y))
grid_mag= np.sqrt(grid_u**2 + grid_v**2 + grid_w**2)

grid_u_inflow=griddata(XY,E1_5_inflow.U,(grid_x,grid_y))
grid_v_inflow=griddata(XY,E1_5_inflow.V,(grid_x,grid_y))
grid_w_inflow=griddata(XY,E1_5_inflow.W,(grid_x,grid_y))
grid_mag_inflow= np.sqrt(grid_u_inflow**2 + grid_v_inflow**2 + grid_w_inflow**2)

grid_u_coarse=griddata(XY,E1_5_coarse.U,(grid_x,grid_y))
grid_v_coarse=griddata(XY,E1_5_coarse.V,(grid_x,grid_y))
grid_w_coarse=griddata(XY,E1_5_coarse.W,(grid_x,grid_y))
grid_mag_coarse=np.sqrt(grid_u_coarse**2 + grid_v_coarse**2 + grid_w_coarse**2)

grid_u_fine=griddata(XY,E1_5_fine.U,(grid_x,grid_y))
grid_v_fine=griddata(XY,E1_5_fine.V,(grid_x,grid_y))
grid_w_fine=griddata(XY,E1_5_fine.W,(grid_x,grid_y))
grid_mag_fine=np.sqrt(grid_u_fine**2 + grid_v_fine**2 + grid_w_fine**2)


#%%
D=0.1
U_ref=6.6
# E1_5.mag = np.sqrt(E1_5.U**2.0  + E1_5.V**2.0)
# E1_5_inflow.mag = np.sqrt(E1_5_inflow.U**2.0  + E1_5_inflow.V**2.0)
# im0 =axes[0].scatter((E1_5.X-0.5)/D,E1_5.Y/D,10, E1_5_inflow.mag/U_ref)
# im1 =axes[1].scatter((E1_5.X-0.5)/D, E1_5.Y/D,10,E1_5.mag/U_ref)

fig, axes = plt.subplots(nrows=1, ncols=2,figsize=(16,6))
im0 =axes[0].contourf((grid_x-0.5)/D,grid_y/D, grid_mag_inflow/U_ref)
axes[0].quiver((grid_x-0.5)/D,grid_y/D,grid_u_inflow,grid_v_inflow, scale=200)

im1 =axes[1].contourf((grid_x-0.5)/D,grid_y/D, grid_mag/U_ref)
axes[1].quiver((grid_x-0.5)/D, grid_y/D,grid_u,grid_v, scale=200)

im0.set_cmap(cmap=plt.get_cmap('jet'))
im1.set_cmap(cmap=plt.get_cmap('jet'))
im0.set_clim(0,1.4)
im1.set_clim(0,1.4)

axes[0].set_title('LES, base inflow')
axes[1].set_title('LES, optimized inflow')
for i in range(0,2):
    axes[i].set(xlim=(-0.05/D,0.15/D), xticks=[-0.5,0,0.5,1,1.5], \
                ylim=(0,0.15/D), yticks=[0,0.5,1,1.5], \
                xlabel='x/D', ylabel='y/D')

fig.colorbar(im0, ax=axes, orientation='vertical')
fig.savefig('../Results/E1_5_inflow_velocity_magnitude.png')

#%%
fig, axes = plt.subplots(nrows=1, ncols=3,figsize=(20,5))

im0 =axes[0].contourf((grid_x-0.5)/D,grid_y/D, grid_mag_coarse/U_ref,levels=np.arange(0,1.4,0.1))
axes[0].quiver((grid_x-0.5)/D,grid_y/D,grid_u_coarse,grid_v_coarse, scale=200)

im1 =axes[1].contourf((grid_x-0.5)/D,grid_y/D, grid_mag/U_ref,levels=np.arange(0,1.4,0.1))
axes[1].quiver((grid_x-0.5)/D, grid_y/D,grid_u,grid_v, scale=200)

im2 =axes[2].contourf((grid_x-0.5)/D,grid_y/D, grid_mag_fine/U_ref,levels=np.arange(0,1.4,0.1))
# im2 =axes[2].scatter(E1_5_coarse.X/D,E1_5_coarse.XD, grid_mag_fine/U_ref)
axes[2].quiver((grid_x-0.5)/D, grid_y/D,grid_u_fine,grid_v_fine, scale=200)

im0.set_cmap(cmap=plt.get_cmap('jet'))
im1.set_cmap(cmap=plt.get_cmap('jet'))
im2.set_cmap(cmap=plt.get_cmap('jet'))
im0.set_clim(0,1.4)
im1.set_clim(0,1.4)
im2.set_clim(0,1.4)

axes[0].set_title('LES, coarse mesh')
axes[1].set_title('LES, base mesh')
axes[2].set_title('LES, fine mesh')
for i in range(0,3):
    axes[i].set(xlim=(-0.05/D,0.15/D), xticks=[-0.5,0,0.5,1,1.5], \
                ylim=(0,0.15/D), yticks=[0,0.5,1,1.5], \
                xlabel='x/D', ylabel='y/D')

fig.colorbar(im0, ax=axes, orientation='vertical')
fig.savefig('../Results/E1_5_grid_velocity_magnitude.png')

#%% each velocity component
fig, axes = plt.subplots(nrows=2, ncols=3,figsize=(15,8))
axes[0][0].scatter((E1_5.X-0.5)/D,E1_5.Y/D,10, E1_5_inflow.U/U_ref)
axes[0][1].scatter((E1_5.X-0.5)/D,E1_5.Y/D,10, E1_5_inflow.V/U_ref)
axes[0][2].scatter((E1_5.X-0.5)/D,E1_5.Y/D,10, E1_5_inflow.W/U_ref)

axes[1][0].scatter((E1_5.X-0.5)/D,E1_5.Y/D,10, E1_5.U/U_ref)
axes[1][1].scatter((E1_5.X-0.5)/D,E1_5.Y/D,10, E1_5.V/U_ref)
axes[1][2].scatter((E1_5.X-0.5)/D,E1_5.Y/D,10, E1_5.W/U_ref)

for i in range(0,3):
    axes[0][i].set(xlim=(-0.05/D,0.15/D), xticks=[-0.5,0,0.5,1,1.5], \
                ylim=(0,0.15/D), yticks=[0,0.5,1,1.5], \
                xlabel='x/D', ylabel='y/D')
    axes[1][i].set(xlim=(-0.05/D,0.15/D), xticks=[-0.5,0,0.5,1,1.5], \
                ylim=(0,0.15/D), yticks=[0,0.5,1,1.5], \
                xlabel='x/D', ylabel='y/D')


#%% each velocity component, grid sensitivity
fig, axes = plt.subplots(3,3,figsize=(16,16))
axes[0][0].contourf((grid_x-0.5)/D,grid_y/D, grid_u_coarse/U_ref)
axes[1][0].contourf((grid_x-0.5)/D,grid_y/D, grid_v_coarse/U_ref)
axes[2][0].contourf((grid_x-0.5)/D,grid_y/D, grid_w_coarse/U_ref)

axes[0][1].contourf((grid_x-0.5)/D,grid_y/D, grid_u/U_ref)
axes[1][1].contourf((grid_x-0.5)/D,grid_y/D, grid_v/U_ref)
axes[2][1].contourf((grid_x-0.5)/D,grid_y/D, grid_w/U_ref)

axes[0][2].contourf((grid_x-0.5)/D,grid_y/D, grid_u_fine/U_ref)
axes[1][2].contourf((grid_x-0.5)/D,grid_y/D, grid_v_fine/U_ref)
axes[2][2].contourf((grid_x-0.5)/D,grid_y/D, grid_w_fine/U_ref)



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



#%%
 # 0.0359375 , 0.0365625 , 0.0371875 ,
       # 0.0378125 , 0.0384375 , 0.0390625 , 0.0396875 , 0.0403125 ,
       # 0.0409375 , 0.0415625 , 0.0421875 
       
y_target = 0.0403125
# y_target = 0.0415625
# y_target = 0.0396875
idx1 = (E1_5_coarse.Y == y_target) 
idx2 = (E1_5.Y == y_target) 

idx3 = (E1_5_fine.Y == y_target) 

plt.plot(data_ref[:,0],data_ref[:,1],'k.',markersize=10)
plt.plot((E1_5_coarse.X[idx1]-0.5)/D, E1_5_coarse.U[idx1]/Uref,'b',linewidth=2)
plt.plot((E1_5.X[idx2]-0.5)/D, E1_5.U[idx2]/Uref,'r',linewidth=2)
plt.plot((E1_5_fine.X[idx3]-0.5)/D, E1_5_fine.U[idx3]/Uref,'g',linewidth=2)

plt.legend(['PIV','LES: coarse mesh','LES: base mesh','LES: fine mesh'])
plt.xlim([-0.25, 1.25])
plt.ylim([0, 1])
plt.xticks([-0.25,0,1,1.25])
plt.xlabel('x/D')
plt.ylabel('U/U$_{ref}$')
plt.grid()
plt.savefig('../Results/E1_5_grid_sensitivity.png')



