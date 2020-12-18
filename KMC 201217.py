# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 20:45:03 2020

@author: Chang
"""

import numpy as np
import matplotlib.pyplot as plt
                                         
Ca0 = 400
Cb0 = 0
Cc0 = 0
CT = Ca0+Cb0+Cc0

k1 = 1.5
k2 = 1 
k = np.array([k1,k2,0])

N = 400
Na = N*Ca0/CT
Nb = N*Cb0/CT
Nc = N*Cc0/CT
Vol = N/CT
                                         # Initial concentration matrix
Nj = np.array([Na,Nb,Nc])                   # for rate calculation
C = np.array([Ca0,Cb0,Cc0])
CC = np.array([Ca0,Cb0,Cc0])                   # for saving data
                                            
t = 0 
T = np.array([0])
                                            # Stoichiometric V 
V = np.array([[-1,0],[1,-1],[0,1]]) 
np.random.seed(260) 

Ca = C[0]
Cb = C[1]
while Ca>=0 and Cb>=0:
    r = k*C
    Rt = r[0]+r[1]
    P1 = np.random.rand()                  # generate random number between(0,1)
    tau = -1*np.log(P1)/Rt          # calculate time interval
    t = t+tau                        # t(1) = t(0)+tau
    P2= np.random.rand()
    a = r[0]/Rt
    
    if P2 < a:
        Nj = Nj+V[:,0]
    else:
        Nj = Nj+V[:,1]
    
    C = Nj/Vol
    CC =np.c_[CC,C]
    Ca = C[0]
    Cb = C[1] 
    T =np.r_[T,[t]]

xx = np.transpose(CC)

plt.plot(T,xx)
plt.xlabel("Time",fontsize=20)
plt.ylabel("C",fontsize=20)
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
plt.xlim(0,10) 
plt.ylim(0,Ca0)
plt.show()

# Analytical solution
tt = np.arange(0,10,0.1)
Ca = Ca0*np.exp(-k1*tt)
Cb = (k1*Ca0/(k2-k1))*(np.exp(-k1*tt)-np.exp(-k2*tt))
Cc = Ca0 - Ca - Cb
plt.plot(tt,Ca)
plt.plot(tt,Cb)
plt.plot(tt,Cc)
plt.show()  