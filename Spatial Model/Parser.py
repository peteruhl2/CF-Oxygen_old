# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 17:21:05 2021

@author: peter
"""

# import matplotlib.pyplot as plt
# import seaborn as sns
import numpy as np
import os
import pandas as pd

# get to folder for a simulation
path = r'C:\Users\peter\OneDrive\Documents\GitHub\CF-Oxygen\Spatial Model\Simulations'
sim = r'Sim g = 0.75, k = 1600'
path += '\\' + sim
os.chdir(path)

print(os.getcwd())

# allocate memory
n_sims = 50
t_steps = 6000 # just needs to be long enough

C = np.zeros((t_steps,n_sims))
F = np.zeros((t_steps,n_sims))
C_avg = np.zeros((t_steps,1))
F_avg = np.zeros((t_steps,1))
C_25 = np.zeros((t_steps,1))
F_25 = np.zeros((t_steps,1))
C_975 = np.zeros((t_steps,1))
F_975 = np.zeros((t_steps,1))

# array for max times
max_times = []

# loop to read in simulation data
for i in range(1,n_sims + 1):
    # sim number
    n = str(i)
    
    # set file name and read in as numpy array
    file_name = 'sim_' + n + '.csv'
    sim = pd.read_csv(file_name, delimiter=',', header=None)
    sim = np.asarray(sim)
    
    # length of current run
    max_t = len(sim[:,0])
    max_times.append(max_t)
    
    # copy to C and F arrays
    C[:max_t,i-1] = sim[:,0]
    F[:max_t,i-1] = sim[:,1]
    
    
# loop to get averages and percentiles
# and write stuff to file
with open('AB_avgs.csv', 'w') as f:
    for i in range(len(C_avg)):
        # means
        # C_avg[i,0] = np.mean(C[i,:])
        # F_avg[i,0] = np.mean(F[i,:])
        C_avg[i,0] = C[i,:][np.nonzero(C[i,:])].mean()
        F_avg[i,0] = F[i,:][np.nonzero(F[i,:])].mean()
        
        # 2.5th percentiles 
        C_25[i,0] = np.percentile(C[i,:][np.nonzero(C[i,:])], q = 2.5)
        F_25[i,0] = np.percentile(F[i,:][np.nonzero(F[i,:])], q = 2.5)
        
        # 97.5th percentiles 
        C_975[i,0] = np.percentile(C[i,:][np.nonzero(C[i,:])], q = 97.5)
        F_975[i,0] = np.percentile(F[i,:][np.nonzero(F[i,:])], q = 97.5)
    
        # comment this if things don't work
        # f.write('%d,%d,%d,%d,%d,%d\n' % (C_25[i],C_avg[i],C_975[i],F_25[i],F_avg[i],F_975[i]))
        f.write('%f,%f,%f,%f,%f,%f\n' % (C_25[i],C_avg[i],C_975[i],F_25[i],F_avg[i],F_975[i]))