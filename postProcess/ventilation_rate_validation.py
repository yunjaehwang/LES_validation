#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 15:29:57 2021

@author: yunjaeh
"""

import numpy as np

data_temp = np.loadtxt('../LES_ventilation/E1_5/output/IN_u.dat')

Time = data_temp[:,0]
Inlet_U = data_temp[:,5]

data_temp = np.loadtxt('../LES_ventilation/E1_5/output/OUT_u.dat')
Outlet_U = data_temp[:,5]

#%%
data_ref = np.loadtxt('../ReferenceData/E1_nondimensional_velocity.csv',\
                      skiprows=2,delimiter=',')

U_PIV = data_ref[0][1]
U_HotFilm = data_ref[0][3]
#%%

dt=0.0001
steps=np.arange(0,40000,10)
dsteps=steps*dt

U_avg=(Inlet_U[1:,]+Outlet_U)/2
U_ref = 6.6

fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(10,4))
# axes[0].plot(Inlet_U/U_ref,'b')
# axes[0].plot(Outlet_U/U_ref,'r')
axes[0].plot(dsteps,U_avg/U_ref,'k')
axes[0].plot([0, 4], [U_PIV, U_PIV],'b--')
axes[0].plot([0, 4], [U_HotFilm, U_HotFilm],'r--')
axes[0].set(xlim=(0,4), xlabel='Time [sec]',\
            ylim=(0,1), ylabel='$U/U_{ref}$')

    
# axes[1].hist(Inlet_U/U_ref,50,alpha=0.5)
# axes[1].hist(Outlet_U/U_ref,50,alpha=0.5)
axes[1].hist(U_avg/U_ref,50,color='gray',alpha=0.7,density=True)
axes[1].plot([U_PIV, U_PIV],[0, 6], 'b--')
axes[1].plot([U_HotFilm, U_HotFilm],[0, 6], 'r--')
axes[1].set(xlim=(0,1), xlabel='$U/U_{ref}$',\
            ylim=(0,6), ylabel='Probability')
fig.legend(['LES','Measurement: PIV', 'Measurement: Hot Film'],\
           bbox_to_anchor=(0.75,0),ncol=3)
fig.savefig('../Results/ventilation_rate_validation.png')





