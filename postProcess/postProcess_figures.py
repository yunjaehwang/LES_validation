#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 20:22:59 2021

@author: yunjaeh
"""

#%% 
import numpy as np
import matplotlib.pyplot as plt

#%%
class fCase():
    def __init__(self):
        self.XYZ =[]
        self.U = []
        self.V = []
        self.W = []
        self.P = []
        self.C = []
        self.C0 = 100
        self.AoA =[]
        self.dt = 0.0001
        self.dstep = 100
        self.fPath =''
    
    def readCoord(self):
        self.XYZ = np.loadtxt(self.fPath+'_XYZ.dat',delimiter=',')
    def readVelocity(self):
        self.U   = np.loadtxt(self.fPath+'_U.dat',delimiter=',')
        self.V   = np.loadtxt(self.fPath+'_V.dat',delimiter=',')
        self.W   = np.loadtxt(self.fPath+'_W.dat',delimiter=',')
    def readP(self):
        self.P   = np.loadtxt(self.fPath+'_P.dat',delimiter=',')
    def readC(self):
        self.C   = np.loadtxt(self.fPath+'_C.dat',delimiter=',')
    def computeAoA(self):
        self.AoA = np.trapz(self.C, dx=self.dt*self.dstep,axis=0)
            

#%%
E1_5, E1_10, E1_20 = fCase(), fCase(), fCase()
E1_5.fPath='../Results/E1_5'
E1_5.readC()
E1_5.computeAoA()

E1_10.fPath='../Results/E1_10'
E1_10.readC()
E1_10.computeAoA()

E1_20.fPath='../Results/E1_20'
E1_20.readC()
E1_20.computeAoA()

E1_5_45=fCase()
E1_5_45.fPath='../Results/E1_5.45'
E1_5_45.readC()
E1_5_45.computeAoA()

E1_5_90=fCase()
E1_5_90.fPath='../Results/E1_5.90'
E1_5_90.readC()
E1_5_90.computeAoA()

#%%
kwargs = dict(histtype='stepfilled', density=True, bins=100, \
              color='b', edgecolor='b', alpha=0.3)

plt.figure(figsize=(6,4))
plt.hist(E1_5.AoA, **kwargs, label='Wall porosity =  5%')
plt.hist(E1_10.AoA,hatch='/', **kwargs, label='Wall porosity = 10%')
plt.hist(E1_20.AoA,hatch='+', **kwargs, label='Wall porosity = 20%')
plt.xlabel('Age of air [second]')
plt.ylabel('Probability')
plt.xlim(0,1)
plt.ylim(0,15)
plt.legend()
plt.savefig('../Results/E1_wall_porosity.png')
    
#%% E1, different angles
kwargs = dict(histtype='stepfilled', density=True, bins=100, \
              edgecolor='k',color='b')

plt.figure(figsize=(6,4))
plt.hist(E1_5.AoA,    alpha=0.3, **kwargs, label='Flow angle =  0$^\circ$')
plt.hist(E1_5_45.AoA, alpha=0.5, **kwargs, label='Flow angle = 45$^\circ$')
plt.hist(E1_5_90.AoA, alpha=0.9, **kwargs, label='Flow angle = 90$^\circ$')
plt.xlabel('Age of air [second]')
plt.ylabel('Probability')
plt.xlim(0,2.5)
plt.ylim(0,10)
plt.legend()
plt.title('Opening locations: center / center')
plt.savefig('../Results/E1_flow_directions.png')

#%%
A1_5=fCase()
A1_5.fPath='../Results/A1_5'
A1_5.readC(),   A1_5.computeAoA()

C1_5=fCase()
C1_5.fPath='../Results/C1_5'
C1_5.readC(),   C1_5.computeAoA()

D1_5=fCase()
D1_5.fPath='../Results/D1_5'
D1_5.readC(),   D1_5.computeAoA()

A1_5_45=fCase()
A1_5_45.fPath='../Results/A1_5.45'
A1_5_45.readC(),   A1_5_45.computeAoA()

A1_5_90=fCase()
A1_5_90.fPath='../Results/A1_5.90'
A1_5_90.readC(),   A1_5_90.computeAoA()

C1_5_45=fCase()
C1_5_45.fPath='../Results/C1_5.45'
C1_5_45.readC(),   C1_5_45.computeAoA()

C1_5_90=fCase()
C1_5_90.fPath='../Results/C1_5.90'
C1_5_90.readC(),   C1_5_90.computeAoA()


D1_5_45=fCase()
D1_5_45.fPath='../Results/D1_5.45'
D1_5_45.readC(),   D1_5_45.computeAoA()

D1_5_90=fCase()
D1_5_90.fPath='../Results/D1_5.90'
D1_5_90.readC(),   D1_5_90.computeAoA()

    
#%%
kwargs = dict(histtype='stepfilled', density=True, bins=40, \
              edgecolor='k',alpha=0.3)

plt.figure(figsize=(6,4))
plt.hist(E1_5.AoA, color='b', **kwargs, label='Opening locations : center / center')
plt.hist(A1_5.AoA, color='r', **kwargs, label='Opening locations : top / top')
plt.hist(C1_5.AoA, color='g', **kwargs, label='Opening locations : bottom / bottom')
plt.hist(D1_5.AoA, color='y', **kwargs, label='Opening locations : top / bottom')
plt.xlabel('Age of air [second]')
plt.ylabel('Probability')
plt.xlim(0,1.)
plt.ylim(0,14)
plt.legend()
# plt.legend( bbox_to_anchor=(0.8,-0.2))
plt.title('Wall porosity = 5%, flow angle = 0$^\circ$')
plt.savefig('../Results/configurations_5%.png')

#%% top / top
kwargs = dict(histtype='stepfilled', density=True, bins=100, \
              edgecolor='k',color='r')

plt.figure(figsize=(6,4))
plt.hist(A1_5.AoA,    alpha=0.3, **kwargs, label='Flow angle =  0$^\circ$')
plt.hist(A1_5_45.AoA, alpha=0.5, **kwargs, label='Flow angle = 45$^\circ$')
plt.hist(A1_5_90.AoA, alpha=0.9, **kwargs, label='Flow angle = 90$^\circ$')
plt.xlabel('Age of air [second]')
plt.ylabel('Probability')
plt.xlim(0,2.5)
plt.ylim(0,10)
plt.legend()
plt.title('Opening locations: top / top')
plt.savefig('../Results/A1_flow_directions.png')

#%% bottom / bottom 
kwargs = dict(histtype='stepfilled', density=True, bins=100, \
              edgecolor='k',color='g')

plt.figure(figsize=(6,4))
plt.hist(C1_5.AoA,    alpha=0.3, **kwargs, label='Flow angle =  0$^\circ$')
plt.hist(C1_5_45.AoA, alpha=0.5, **kwargs, label='Flow angle = 45$^\circ$')
plt.hist(C1_5_90.AoA, alpha=0.9, **kwargs, label='Flow angle = 90$^\circ$')
plt.xlabel('Age of air [second]')
plt.ylabel('Probability')
plt.xlim(0,2.5)
plt.ylim(0,10)
plt.legend()
plt.title('Opening locations: bottom / bottom')
plt.savefig('../Results/C1_flow_directions.png')

#%% top / bottom

kwargs = dict(histtype='stepfilled', density=True, bins=100, \
              edgecolor='k',color='y')

plt.figure(figsize=(6,4))
plt.hist(D1_5.AoA,    alpha=0.3, **kwargs, label='Flow angle =  0$^\circ$')
plt.hist(D1_5_45.AoA, alpha=0.5, **kwargs, label='Flow angle = 45$^\circ$')
plt.hist(D1_5_90.AoA, alpha=0.9, **kwargs, label='Flow angle = 90$^\circ$')
plt.xlabel('Age of air [second]')
plt.ylabel('Probability')
plt.xlim(0,2.5)
plt.ylim(0,10)
plt.legend()
plt.title('Opening locations: top / bottom')
plt.savefig('../Results/D1_flow_directions.png')

#%% 45 cases
kwargs = dict(histtype='stepfilled', density=True, bins=100, \
              edgecolor='k',alpha=0.5)

plt.figure(figsize=(6,4))
plt.hist(E1_5_45.AoA, color='b', **kwargs, label='Opening locations : center / center')
plt.hist(A1_5_45.AoA, color='r', **kwargs, label='Opening locations : top / top')
plt.hist(C1_5_45.AoA, color='g', **kwargs, label='Opening locations : bottom / bottom')
plt.hist(D1_5_45.AoA, color='y', **kwargs, label='Opening locations : top / bottom')
plt.xlabel('Age of air [second]')
plt.ylabel('Probability')
plt.xlim(0,1.2)
plt.ylim(0,12)
plt.legend()
# plt.legend( bbox_to_anchor=(0.8,-0.2))
plt.title('Wall porosity = 5%, flow angle = 45$^\circ$')
plt.savefig('../Results/configurations_45deg.png')

#%% 90 degree
kwargs = dict(histtype='stepfilled', density=True, bins=100, \
              edgecolor='k',alpha=0.7)

plt.figure(figsize=(6,4))
plt.hist(E1_5_90.AoA, color='b', **kwargs, label='Opening locations : center / center')
plt.hist(A1_5_90.AoA, color='r', **kwargs, label='Opening locations : top / top')
plt.hist(C1_5_90.AoA, color='g', **kwargs, label='Opening locations : bottom / bottom')
plt.hist(D1_5_90.AoA, color='y', **kwargs, label='Opening locations : top / bottom')
plt.xlabel('Age of air [second]')
plt.ylabel('Probability')
plt.xlim(0,2.5)
plt.ylim(0,6)
plt.legend()
plt.title('Wall porosity = 5%, flow angle = 90$^\circ$')
plt.savefig('../Results/configurations_90deg.png')

#%% print data

print('Mean')
print('E1', np.mean(E1_5.AoA), np.mean(E1_5_45.AoA), np.mean(E1_5_90.AoA))
print('A1', np.mean(A1_5.AoA), np.mean(A1_5_45.AoA), np.mean(A1_5_90.AoA))
print('C1', np.mean(C1_5.AoA), np.mean(C1_5_45.AoA), np.mean(C1_5_90.AoA))
print('D1', np.mean(D1_5.AoA), np.mean(D1_5_45.AoA), np.mean(D1_5_90.AoA))

print('Std')
print('E1', np.std(E1_5.AoA), np.std(E1_5_45.AoA), np.std(E1_5_90.AoA))
print('A1', np.std(A1_5.AoA), np.std(A1_5_45.AoA), np.std(A1_5_90.AoA))
print('C1', np.std(C1_5.AoA), np.std(C1_5_45.AoA), np.std(C1_5_90.AoA))
print('D1', np.std(D1_5.AoA), np.std(D1_5_45.AoA), np.std(D1_5_90.AoA))






















    



