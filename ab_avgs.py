# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 08:55:47 2020

@author: peter
"""

# import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
# import seaborn as sns

path = r'C:\Users\peter\OneDrive\Desktop\cyst fib\OxygenModels\Simulations'
os.chdir(path)

n_sims = 2000
t_steps = 960

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
    
    # copy to C and F arrays
    C[:,i-1] = sim[:,0]
    F[:,i-1] = sim[:,1]
    

# loop to get averages and percentiles
# and write stuff to file
with open('AB_avgs.csv', 'w') as f:
    for i in range(len(C_avg)):
        # means
        C_avg[i,0] = np.mean(C[i,:])
        F_avg[i,0] = np.mean(F[i,:])
        
        # 2.5th percentiles 
        C_25[i,0] = np.percentile(C[i,:],q = 2.5)
        F_25[i,0] = np.percentile(F[i,:],q = 2.5)
        
        # 97.5th percentiles 
        C_975[i,0] = np.percentile(C[i,:],q = 97.5)
        F_975[i,0] = np.percentile(F[i,:],q = 97.5)
    
        # comment this if things don't work
        f.write('%d,%d,%d,%d,%d,%d\n' % (C_25[i],C_avg[i],C_975[i],F_25[i],F_avg[i],F_975[i]))
    


    