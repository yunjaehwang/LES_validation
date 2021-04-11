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
fig, axes = plt.subplots(nrows=1, ncols=2,figsize=(10,6))
axes[0].scatter(E1_5.X,E1_5.Y,10,E1_5.U)
axes[1].scatter(E1_5_inflow.X, E1_5_inflow.Y,10,E1_5_inflow.U)

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

plt.savefig('../Results/')




