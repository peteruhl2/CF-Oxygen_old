# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 08:55:47 2020
Reused June 2021

@author: peter

comment short cut: Ctrl + 1
"""

# import matplotlib.pyplot as plt
# import seaborn as sns
import numpy as np
import os
import pandas as pd

# path = r'C:\Users\peter\OneDrive\Desktop\cyst fib\OxygenModels\Simulations' # from Dec 2021
path = r'C:\Users\peter\OneDrive\Documents\GitHub\CF-Oxygen\Non-spatial Model\Simulations'
os.chdir(path)

n_sims = 500
t_steps = 1440

C = np.zeros((t_steps,n_sims))
F = np.zeros((t_steps,n_sims))
C_avg = np.zeros((t_steps,1))
F_avg = np.zeros((t_steps,1))
C_25 = np.zeros((t_steps,1))
F_25 = np.zeros((t_steps,1))
C_975 = np.zeros((t_steps,1))
F_975 = np.zeros((t_steps,1))

# loop to read in simulation data
for i in range(1,n_sims + 1):
    # sim number
    n = str(i)
    
    # set file name and read in as numpy array
    file_name = 'sim_' + n + '.csv'
    sim = pd.read_csv(file_name, delimiter=' ', header=None)
    sim = np.asarray(sim)
    
    # # copy to C and F arrays
    # C[:,i-1] = sim[:,0]
    # F[:,i-1] = sim[:,1]
    
    # copy to C and F arrays
    C[:,i-1] = sim[:,0]/sim[:,3]
    F[:,i-1] = sim[:,1]/sim[:,3]
    
    # replace 0's in F with nan
    C[F == 0] = np.nan
    F[F == 0] = np.nan
    

# loop to get averages and percentiles
# and write stuff to file
with open('AB_avgs.csv', 'w') as f:
    for i in range(len(C_avg)):
        # means
        C_avg[i,0] = np.nanmean(C[i,:])
        F_avg[i,0] = np.nanmean(F[i,:])
        
        # 2.5th percentiles 
        C_25[i,0] = np.nanpercentile(C[i,:],q = 2.5)
        F_25[i,0] = np.nanpercentile(F[i,:],q = 2.5)
        
        # 97.5th percentiles 
        C_975[i,0] = np.nanpercentile(C[i,:],q = 97.5)
        F_975[i,0] = np.nanpercentile(F[i,:],q = 97.5)
    
        # comment this if things don't work
        f.write('%d,%d,%d,%d,%d,%d\n' % (C_25[i],C_avg[i],C_975[i],F_25[i],F_avg[i],F_975[i]))
            
