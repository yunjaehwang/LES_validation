#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 14:41:50 2021

@author: yunjaeh
"""

import numpy as np
import matplotlib.pyplot as plt

fCase = 'D1_5.90'
# fPath = '../LES_ventilation/E1_5.base_inflow/output/PC/'
fPath = '../LES_ventilation/'+fCase+'/post/PC/'

coord = np.loadtxt(fPath+'C.pxyz', skiprows=1)


dt, dstep = 0.0001, 100

C0 = 100
rawC, rawU, rawV, rawW, rawP = [], [], [], [], []

steps     = np.arange(0,90001,dstep)
timeSteps = steps*dt
timeSteps = timeSteps-timeSteps[0]


for i in steps:
    print('Reading: ', fPath+'C.'+str(i).zfill(8)+'.pcd')
    data=np.loadtxt(fPath+'C.'+str(i).zfill(8)+'.pcd',skiprows=1)
    rawC.append(data[:,0])
    rawU.append(data[:,1])
    rawV.append(data[:,2])
    rawW.append(data[:,3])
    rawP.append(data[:,4])

# choose only points inside the house and save data
idx = np.asarray(rawC)[0,:] > 90.0
XYZ = coord[idx,1:]
C = np.asarray(rawC)[:,idx]/C0
U = np.asarray(rawU)[:,idx]
V = np.asarray(rawV)[:,idx]
W = np.asarray(rawW)[:,idx]
P = np.asarray(rawP)[:,idx]
Umag = np.sqrt(U**2+V**2+W**2)

[lenC, lenPt] = np.shape(C)

np.savetxt('../Results/'+fCase+'_XYZ.dat',XYZ, delimiter=',', fmt='%4e')
np.savetxt('../Results/'+fCase+'_U.dat',U, delimiter=',', fmt='%4e')
np.savetxt('../Results/'+fCase+'_V.dat',V, delimiter=',', fmt='%4e')
np.savetxt('../Results/'+fCase+'_W.dat',W, delimiter=',', fmt='%4e')
np.savetxt('../Results/'+fCase+'_C.dat',C, delimiter=',', fmt='%4e')
np.savetxt('../Results/'+fCase+'_P.dat',P, delimiter=',', fmt='%4e')


#%%
C_mean = np.mean(C,axis=0)
plt.figure()
plt.plot(C_mean,XYZ[:,1],'k.')
# plt.xlim(21,31)
# plt.xlim(16,21)
# axes[0].set(xlabel='Time [sec]', ylabel='C/C0 [-]', \
            # title='Time evolution of concentration')

#%%
idY0 = (XYZ[:,1] == 0.06)
idZ0 = (XYZ[:,2] == 0.55)

# idY1 = (XYZ_e1_5_opt[:,1] == 0.04)
# idZ1 = (XYZ_e1_5_opt[:,2] == 0.55)


plt.scatter(XYZ[idZ0,0], XYZ[idZ0,1], c=np.mean(U,axis=0)[idZ0])

# fig, axes = plt.subplots(nrows=1, ncols=3)
# axes[0].scatter(XYZ[idZ0,], XYZ[idZ0])
# axes[1].scatter()
# axes[2].scatter()

#%%
plt.plot(XYZ[idY0&idZ0,0], np.mean(U,axis=0)[idY0&idZ0]/6.6)
# plt.plot(XYZ_e1_5_opt[idY1&idZ1,0], np.mean(U_e1_5_opt,axis=0)[idY1&idZ1]/6.6)
# plt.title('E1_5, base_inflow')



#%% age of air calculation and plot
# aoa=[]
# for i in range(1,lenC):
    # aoa.append(np.trapz(C[0:i,:],dx=dt*dstep,axis=0))
# aoa=np.asarray(aoa)

aoa = np.trapz(C[0:i,:],dx=dt*dstep,axis=0)
   
# #%%
# plt.figure()
# for i in range(0,lenPt,100):
#     plt.plot(timeSteps[1:], aoa[:,i],'gray')
# plt.xlabel('Time [sec]')
# plt.ylabel('Age of air [sec]')
# plt.title('Time evolution of age of air at indoor locations')

#%%
plt.scatter(XYZ[idY0,0], XYZ[idY0,2],c=aoa[idY0])
plt.colorbar()

#%%
plt.figure()
# plt.hist(aoa,25, density=True)
plt.hist(aoa,100, density=True)
# plt.hist(aoa[-10,:],50, density=True)
# plt.hist(aoa[-100,:],50, density=True)
plt.xlabel('Age of air [sec]')
plt.ylabel('Density')
plt.title('PDF of age of air')

#%% comparison between two1


# %%

# plt.figure()
plt.hist(3600/aoa,density=True,bins=range(0,10000,100))
# plt.xlabel('ACH [1/hr]')
# plt.ylabel('Density')
# plt.title('PDF of age of air')
# plt.xlim(0,40)

#%%
# plt.figure()
# plt.hist(aoa_E1_5,25,facecolor='blue',density=True, alpha=0.5)
# plt.hist(aoa_A1_5,25,facecolor='red',density=True, alpha=0.5)
# plt.hist(aoa_C1_5,25,facecolor='green',density=True, alpha=0.5)
# plt.hist(aoa_D1_5,25,facecolor='yellow',density=True, alpha=0.5)



#%% remaing Concentration
fig, axes = plt.subplots(ncols=2, nrows=1)

plt.figure()
for i in range(0,lenPt,100):
    axes[0].plot(timeSteps, C[:,i],'gray')
axes[0].plot(timeSteps, np.mean(C,axis=1),'k')

axes[0].set(xlabel='Time [sec]', ylabel='C/C0 [-]', \
            title='Time evolution of concentration')

axes[1].hist(C[-1,:])


#%% read image
data_fine    =np.loadtxt('../LES_ventilation/E1_5.fine/output/images_stats_old/mag.avg.raw_values.dat',skiprows=1);
data_fine_new=np.loadtxt('../LES_ventilation/E1_5.fine/output/images_stats_new/mag.avg.raw_values.dat',skiprows=1);
data_fine_all=np.loadtxt('../LES_ventilation/E1_5.fine/output/images_stats_all/mag.avg.raw_values.dat',skiprows=1);
data_base    =np.loadtxt('../LES_ventilation/E1_5/output/images_stats/mag.avg.raw_values.dat',skiprows=1);
data_inflow  =np.loadtxt('../LES_ventilation/E1_5.base_inflow/output/images_stats/mag.avg.raw_values.dat',skiprows=1);    

# %%
# plt.scatter(data_png[:,4], data_png[:,5], c=data_png[:,7])
idY_base = (data_base[:,5] < 0.04) & (data_base[:,5] > 0.039)
idY_fine = (data_fine[:,5] < 0.04) & (data_fine[:,5] > 0.039)
idY_fine_new = (data_fine_new[:,5] < 0.04) & (data_fine_new[:,5] > 0.039)
# idY_fine_all = (data_fine_all[:,5] < 0.04) & (data_fine_all[:,5] > 0.039)
# idY_inflow = (data_inflow[:,5] < 0.041) & (data_inflow[:,5] > 0.039)

plt.plot(data_base[idY_base,4], data_base[idY_base,7]/6.6,'k')
plt.plot(data_fine[idY_fine,4], data_fine[idY_fine,7]/6.6,'r')
plt.plot(data_fine_new[idY_fine_new,4], data_fine_new[idY_fine_new,7]/6.6,'g')
# plt.plot(data_fine_all[idY_fine_all,4], data_fine_all[idY_fine_all,7]/6.6,'b')
# plt.plot(data_inflow[idY_inflow,4], data_inflow[idY_inflow,7]/6.6,'b')
plt.xlim(0.475, 0.625)
# plt.plot(XYZ[idY0&idZ0,0], np.mean(Umag,axis=0)[idY0&idZ0]/6.6)







