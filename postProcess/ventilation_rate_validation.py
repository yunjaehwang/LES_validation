#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 11 15:29:57 2021

@author: yunjaeh
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class fCase():
    def __init__(self, amplification_factor=1.0, wall_porosity=5):
        self.time=[]
        self.step=[]
        self.U_inlet=[]
        self.U_outlet=[]
        self.U_avg=[]
        self.U_reshaped = None
        self.P_inlet =[]
        self.P_outlet =[]
        self.deltaP =[]
        self.dt = 0.0001
        self.dstep = 10
        
        self.rho = 1.225
        self.mu = 1.5*10**-5
        h, w = 0.018, 0.023*wall_porosity/5
        self.d_h = 4*w*h/(2*w+2*h)
        
        self.amplification_factor = amplification_factor
        self.wall_porosity = wall_porosity
    def readData(self,fPath):
        data_temp = np.loadtxt(fPath+'IN_u.dat')
        self.step = data_temp[:,0]
        self.time = self.step*self.dt*self.dstep
        
        self.U_inlet = data_temp[:,5] * self.amplification_factor
        data_temp = np.loadtxt(fPath+'OUT_u.dat')
        self.U_outlet = data_temp[:,5] * self.amplification_factor
        self.U_avg = (self.U_inlet + self.U_outlet)/2.0

        data_temp = np.loadtxt(fPath+'IN_p.dat')
        self.P_inlet = data_temp[:,5]
        data_temp = np.loadtxt(fPath+'OUT_p.dat')
        self.P_outlet = data_temp[:,5]
        self.deltaP = np.abs(self.P_inlet - self.P_outlet) 
    def reshape(self,length):
        self.U_reshaped = np.mean(self.U_avg[1:].reshape((-1,length)),axis=1)
        self.deltaP_reshaped = np.mean(self.deltaP[1:].reshape((-1,length)),axis=1)
    def computeQTY(self):
        self.Cd =  self.U_avg / np.sqrt(2*np.abs(self.deltaP)/self.rho)
        self.Re_o = self.U_avg * self.d_h *self.rho / self.mu
        if(self.U_reshaped is not None):
            print('Process reshaped data')
            self.Cd_reshaped = self.U_reshaped/np.sqrt(2*np.abs(self.deltaP_reshaped)/self.rho) 
            self.Re_o_reshaped = self.U_reshaped *self.d_h *self.rho /self.mu



#%%
E1_5=fCase()
E1_5.readData('../LES_ventilation/E1_5/post/')
E1_5.computeQTY()

#%% Read reference data
data_ref = np.loadtxt('../ReferenceData/E1_nondimensional_velocity.csv',\
                      skiprows=2,delimiter=',')
U_PIV = data_ref[0][1]
U_HotFilm = data_ref[0][3]
#%%

A1_5, A1_5_45, A1_5_90 = fCase(), fCase(amplification_factor=2), fCase()
C1_5, C1_5_45, C1_5_90 = fCase(), fCase(amplification_factor=2), fCase()
D1_5, D1_5_45, D1_5_90 = fCase(), fCase(amplification_factor=2), fCase()
E1_5, E1_5_45, E1_5_90 = fCase(), fCase(amplification_factor=2), fCase()
E1_10, E1_20 = fCase(wall_porosity=10), fCase(wall_porosity=20)

A1_5.readData('../LES_ventilation/A1_5/post/')
A1_5_45.readData('../LES_ventilation/A1_5.45/post/')
A1_5_90.readData('../LES_ventilation/A1_5.90/post/')

C1_5.readData('../LES_ventilation/C1_5/post/')
C1_5_45.readData('../LES_ventilation/C1_5.45/post/')
C1_5_90.readData('../LES_ventilation/C1_5.90/post/')

D1_5.readData('../LES_ventilation/D1_5/post/')
D1_5_45.readData('../LES_ventilation/D1_5.45/post/')
D1_5_90.readData('../LES_ventilation/D1_5.90/post/')

E1_5.readData('../LES_ventilation/E1_5/post/')
E1_5_45.readData('../LES_ventilation/E1_5.45/post/')
E1_5_90.readData('../LES_ventilation/E1_5.90/post/')

E1_10.readData('../LES_ventilation/E1_10/post/')
E1_20.readData('../LES_ventilation/E1_20/post/')

#%%
U_ref = 6.6

fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(10,4))
# axes[0].plot(Inlet_U/U_ref,'b')
# axes[0].plot(Outlet_U/U_ref,'r')
axes[0].plot(E1_5.time,E1_5.U_avg/U_ref,'k')
axes[0].plot([0, 4], [U_PIV, U_PIV],'b--')
axes[0].plot([0, 4], [U_HotFilm, U_HotFilm],'r--')
axes[0].set(xlim=(0,4), xlabel='Time [sec]',\
            ylim=(0,1), ylabel='$U/U_{ref}$')

    
# axes[1].hist(Inlet_U/U_ref,50,alpha=0.5)
# axes[1].hist(Outlet_U/U_ref,50,alpha=0.5)
axes[1].hist(E1_5.U_avg/U_ref,50,color='gray',alpha=0.7,density=True)
axes[1].plot([U_PIV, U_PIV],[0, 6], 'b--')
axes[1].plot([U_HotFilm, U_HotFilm],[0, 6], 'r--')
axes[1].set(xlim=(0,1), xlabel='$U/U_{ref}$',\
            ylim=(0,6), ylabel='Probability')
fig.legend(['LES','Measurement: PIV', 'Measurement: Hot Film'],\
           bbox_to_anchor=(0.75,0),ncol=3)
# fig.savefig('../Results/ventilation_rate_validation.png')


#%% 
# Q/A = Cd * sqrt(2 delta P / rho)
Cd_ref = 0.61
rho = 1.225

U_p = Cd_ref * np.sqrt(2*np.abs(E1_5.deltaP)/rho)

plt.hist(U_p/U_ref,50, density=True)
plt.hist(E1_5.U_avg/U_ref,50,color='gray',alpha=0.7,density=True)
plt.plot([U_PIV, U_PIV],[0, 6], 'b--')
plt.plot([U_HotFilm, U_HotFilm],[0, 6], 'r--')
# plt.set(xlim=(0,1), xlabel='$U/U_{ref}$',\
            # ylim=(0,6), ylabel='Probability')

#%% backtrack Cd
# Cd = Q / A / sqrt(2 delta P / rho)
#    = U / sqrt(2 delta P / rho)
nAvg = 10
A1_5.reshape(nAvg), A1_5_45.reshape(nAvg),  A1_5_90.reshape(nAvg)
C1_5.reshape(nAvg), C1_5_45.reshape(nAvg),  C1_5_90.reshape(nAvg)
D1_5.reshape(nAvg), D1_5_45.reshape(nAvg),  D1_5_90.reshape(nAvg)
E1_5.reshape(nAvg), E1_5_45.reshape(nAvg),  E1_5_90.reshape(nAvg)
E1_10.reshape(nAvg),E1_20.reshape(nAvg)

A1_5.computeQTY(),  A1_5_45.computeQTY(),   A1_5_90.computeQTY()
C1_5.computeQTY(),  C1_5_45.computeQTY(),   C1_5_90.computeQTY()
D1_5.computeQTY(),  D1_5_45.computeQTY(),   D1_5_90.computeQTY()
E1_5.computeQTY(),  E1_5_45.computeQTY(),   E1_5_90.computeQTY()
E1_10.computeQTY(), E1_20.computeQTY()


kwargs = dict(histtype='stepfilled',density=True, bins=np.arange(0,2.0,0.02), \
           edgecolor='k',alpha=0.5)

fig, ax = plt.subplots(1,1)
A1_hist=ax.hist(A1_5.Cd,color='r',**kwargs)
C1_hist=ax.hist(C1_5.Cd,color='g',**kwargs)
D1_hist=ax.hist(D1_5.Cd,color='y',**kwargs)
E1_hist=ax.hist(E1_5.Cd,color='b',**kwargs)
# plt.hist(E1_10.Cd,bins=hist_bins,density=True)
# plt.hist(E1_20.Cd,bins=hist_bins,density=True)
ax.set(xlim=(0,2.0), xlabel='$C_d$, discharge coefficient [-]', \
       ylabel='Probability');
ax.legend(['Top/Top','Bottom/Bottom','Top/Bottom','Center/Center'])
fig.savefig('../Results/Cd_0deg.png')

fig, ax = plt.subplots(1,1)
A1_45_hist=ax.hist(A1_5_45.Cd,color='r',**kwargs)
C1_45_hist=ax.hist(C1_5_45.Cd,color='g',**kwargs)
D1_45_hist=ax.hist(D1_5_45.Cd,color='y',**kwargs)
E1_45_hist=ax.hist(E1_5_45.Cd,color='b',**kwargs)
# plt.hist(E1_10.Cd,bins=hist_bins,density=True)
# plt.hist(E1_20.Cd,bins=hist_bins,density=True)
ax.set(xlim=(0,2.0), xlabel='$C_d$, discharge coefficient [-]', \
       ylabel='Probability');
ax.legend(['Top/Top','Bottom/Bottom','Top/Bottom','Center/Center'])
fig.savefig('../Results/Cd_45deg.png')

fig, ax = plt.subplots(1,1)
A1_90_hist=ax.hist(A1_5_90.Cd,color='r',**kwargs)
C1_90_hist=ax.hist(C1_5_90.Cd,color='g',**kwargs)
D1_90_hist=ax.hist(D1_5_90.Cd,color='y',**kwargs)
E1_90_hist=ax.hist(E1_5_90.Cd,color='b',**kwargs)
# plt.hist(E1_10.Cd,bins=hist_bins,density=True)
# plt.hist(E1_20.Cd,bins=hist_bins,density=True)
ax.set(xlim=(0,2.0), xlabel='$C_d$, discharge coefficient [-]', \
       ylabel='Probability');
ax.legend(['Top/Top','Bottom/Bottom','Top/Bottom','Center/Center'])
fig.savefig('../Results/Cd_90deg.png')


fig, ax = plt.subplots(1,1)
E1_hist=ax.hist(E1_5.Cd,color='r',**kwargs)
E1_10_hist=ax.hist(E1_10.Cd,color='r',**kwargs)
E1_20_hist=ax.hist(E1_20.Cd,color='r',**kwargs)
ax.set(xlim=(0,2.0), xlabel='$C_d$, discharge coefficient [-]', \
       ylabel='Probability');
fig.savefig('../Results/E1_wallporosity.png')

#%%


kwargs = dict(histtype='stepfilled',density=True, bins=np.arange(0,2.0,0.01), \
           edgecolor='k',alpha=0.5)

fig, axes = plt.subplots(3,4,figsize=(12,7))
axes[0][0].hist(A1_5.Cd,color='b',**kwargs)
axes[0][1].hist(C1_5.Cd,color='g',**kwargs)
axes[0][2].hist(D1_5.Cd,color='y',**kwargs)
axes[0][3].hist(E1_5.Cd,color='r',**kwargs)

axes[1][0].hist(A1_5_45.Cd,color='b',**kwargs)
axes[1][1].hist(C1_5_45.Cd,color='g',**kwargs)
axes[1][2].hist(D1_5_45.Cd,color='y',**kwargs)
axes[1][3].hist(E1_5_45.Cd,color='r',**kwargs)

axes[2][0].hist(A1_5_90.Cd,color='b',**kwargs)
axes[2][1].hist(C1_5_90.Cd,color='g',**kwargs)
axes[2][2].hist(D1_5_90.Cd,color='y',**kwargs)
axes[2][3].hist(E1_5_90.Cd,color='r',**kwargs)

for i in range(0,4):
    axes[0][i].set(ylim=(0,8),yticks=[0,2,4,6,8])
    axes[1][i].set(ylim=(0,8),yticks=[0,2,4,6,8])
    axes[2][i].set(ylim=(0,3),yticks=[0,1,2,3])

axes[0][0].set_title('Top/Top')
axes[0][1].set_title('Bottom/Bottom')
axes[0][2].set_title('Bottom/Top')
axes[0][3].set_title('Center/Center')

axes[0][0].set(ylabel='Direction: $0^\circ$')
axes[1][0].set(ylabel='Direction: $45^\circ$')
axes[2][0].set(ylabel='Direction: $90^\circ$')

fig.savefig('../Results/All_Cd.png')



#%%
bb = []
data_bp={'A1_5':A1_5.Cd, 'A1_5.45':A1_5_45.Cd}

kwargs=dict(width=0.5)
fig, ax = plt.subplots(1,1,figsize=(8,6))
ax.plot([-3, 23],[0.61, 0.61],'k--')
# sns.boxplot(data=[A1_5.Cd, bb, bb, bb, bb, bb, bb, bb, bb,\
#                   A1_5_45.Cd, bb, bb, bb, bb, bb, bb, bb, bb,\
#                   A1_5_90.Cd],orient='v',color='blue')
# sns.boxplot(data=[bb,C1_5.Cd, bb, bb, bb, bb, bb, bb, bb, bb,\
#                   C1_5_45.Cd, bb, bb, bb, bb, bb, bb, bb, bb,\
#                   C1_5_90.Cd],orient='v',color='green')    
# sns.boxplot(data=[bb,bb,D1_5.Cd, bb, bb, bb, bb, bb, bb, bb, bb,\
#                   D1_5_45.Cd, bb, bb, bb, bb, bb, bb, bb, bb,\
#                   D1_5_90.Cd],orient='v',color='y')
# sns.boxplot(data=[bb,bb,bb,E1_5.Cd, bb, bb, bb, bb, bb, bb, bb, bb,\
#                   E1_5_45.Cd, bb, bb, bb, bb, bb, bb, bb, bb,\
#                   E1_5_90.Cd],orient='v',color='red')

# ax.plot([0, 9, 18], [np.median(A1_5.Cd), np.median(A1_5_45.Cd), np.median(A1_5_90.Cd)] ,'b.-', mfc='none')
# ax.plot([1, 10, 19], [np.median(C1_5.Cd), np.median(C1_5_45.Cd), np.median(C1_5_90.Cd)] ,'g.-', mfc='none')
# ax.plot([2, 11, 20], [np.median(D1_5.Cd), np.median(D1_5_45.Cd), np.median(D1_5_90.Cd)] ,'y.-', mfc='none')
# ax.plot([3, 12, 21], [np.median(E1_5.Cd), np.median(E1_5_45.Cd), np.median(E1_5_90.Cd)] ,'r.-', mfc='none')


ax = sns.boxplot(data=[A1_5.Cd_reshaped, bb, bb, bb, bb, bb, bb, bb, bb,\
                  A1_5_45.Cd_reshaped, bb, bb, bb, bb, bb, bb, bb, bb,\
                  A1_5_90.Cd_reshaped],orient='v',color='b',**kwargs)
ax = sns.boxplot(data=[bb,C1_5.Cd_reshaped, bb, bb, bb, bb, bb, bb, bb, bb,\
                  C1_5_45.Cd_reshaped, bb, bb, bb, bb, bb, bb, bb, bb,\
                  C1_5_90.Cd_reshaped],orient='v',color='g',**kwargs)    
ax = sns.boxplot(data=[bb,bb,D1_5.Cd_reshaped, bb, bb, bb, bb, bb, bb, bb, bb,\
                  D1_5_45.Cd_reshaped, bb, bb, bb, bb, bb, bb, bb, bb,\
                  D1_5_90.Cd_reshaped],orient='v',color='y',**kwargs)
ax = sns.boxplot(data=[bb,bb,bb,E1_5.Cd_reshaped, bb, bb, bb, bb, bb, bb, bb, bb,\
                  E1_5_45.Cd_reshaped, bb, bb, bb, bb, bb, bb, bb, bb,\
                  E1_5_90.Cd_reshaped],orient='v',color='r',**kwargs)
    

ftsize=15
plt.xticks([1.5, 10.5, 19.5],['0','45','90'],size=ftsize)
plt.yticks([0,0.5,0.61,1.0,1.5],size=ftsize)

plt.xlim(-2,23)
plt.ylim(0,1.5)

plt.xlabel('Flow angle [degree]',size=ftsize)
plt.ylabel('Discharge coefficient, $C_d$ [-]',size=ftsize)

fig.savefig('../Results/Cd_summary.png')

# plt.legend(['Top/Top','Bottom/Bottom','Bottom/Top','Center/Center'])



# %%
hist_bin=np.arange(0,2.0,0.02) 
A1_peak = hist_bin[np.argmax(A1_hist[0])]
C1_peak = hist_bin[np.argmax(C1_hist[0])]
D1_peak = hist_bin[np.argmax(D1_hist[0])]
E1_peak = hist_bin[np.argmax(E1_hist[0])]

A1_45_peak = hist_bin[np.argmax(A1_45_hist[0])]
C1_45_peak = hist_bin[np.argmax(C1_45_hist[0])]
D1_45_peak = hist_bin[np.argmax(D1_45_hist[0])]
E1_45_peak = hist_bin[np.argmax(E1_45_hist[0])]

A1_90_peak = hist_bin[np.argmax(A1_90_hist[0])]
C1_90_peak = hist_bin[np.argmax(C1_90_hist[0])]
D1_90_peak = hist_bin[np.argmax(D1_90_hist[0])]
E1_90_peak = hist_bin[np.argmax(E1_90_hist[0])]

E1_10_peak = hist_bin[np.argmax(E1_10_hist[0])]
E1_20_peak = hist_bin[np.argmax(E1_20_hist[0])]

print(A1_peak, C1_peak, D1_peak, E1_peak)

#%%

plt.plot(np.mean(A1_5.Re_o), A1_peak, 'bo', mfc='none')    
plt.plot(np.mean(C1_5.Re_o), C1_peak, 'go', mfc='none')
plt.plot(np.mean(D1_5.Re_o), D1_peak, 'yo', mfc='none')
plt.plot(np.mean(E1_5.Re_o), E1_peak, 'ro', mfc='none')

plt.plot(np.mean(A1_5_45.Re_o), A1_45_peak, 'bs', mfc='none')    
plt.plot(np.mean(C1_5_45.Re_o), C1_45_peak, 'gs', mfc='none')
plt.plot(np.mean(D1_5_45.Re_o), D1_45_peak, 'ys', mfc='none')
plt.plot(np.mean(E1_5_45.Re_o), E1_45_peak, 'rs', mfc='none')

plt.plot(np.mean(A1_5_90.Re_o), A1_90_peak, 'bx', mfc='none')    
plt.plot(np.mean(C1_5_90.Re_o), C1_90_peak, 'gx', mfc='none')
plt.plot(np.mean(D1_5_90.Re_o), D1_90_peak, 'yx', mfc='none')
plt.plot(np.mean(E1_5_90.Re_o), E1_90_peak, 'rx', mfc='none')

plt.plot(np.mean(E1_10.Re_o), E1_10_peak,'ro', mfc='g')
plt.plot(np.mean(E1_20.Re_o), E1_20_peak,'ro', mfc='r')


plt.xticks(range(0,15001,2500))
plt.xlim(0,10000)
plt.ylim(0,1.6)
plt.xlabel('$Re_o$')
plt.ylabel('$C_d$')

plt.savefig('../Results/Avg_Re_Cd.png')

# plt.legend(['Wall porosity:  5%','Wall porosity: 10%','Wall porosity: 20%'])
# plt.savefig('../Results/Cd_Re_o_E1_porosity.png')

# plt.savefig('../Results/Cd_Re_o_E1_0deg.png')



