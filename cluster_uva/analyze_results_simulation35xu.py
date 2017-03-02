# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 10:18:48 2017

@author: misWin

analyze simulation35xu data
BPSK, TR and GN models, 
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from milp2_xu import *

#%% read data
iteration_history_tr_list = []
iteration_history_gn_list = []
for i in range(30):
    print('Reading results of simulation {}'.format(i))
    file_name = 'results/results35xu/simulation35xu_{}.pkl'.format(i)
    tmp1, tmp2 = read_data(file_name)
    iteration_history_tr_list.append(tmp1)
    iteration_history_gn_list.append(tmp2)
    
#%% extract c and total
ctr = []
cgn = []
ttr = []
tgn = []
for i in range(30):
    print('Extracting data from simulation {}'.format(i))
    ctr.append(extract_history(iteration_history_tr_list[i], 'c'))
    cgn.append(extract_history(iteration_history_gn_list[i], 'c'))
    ttr.append(extract_history(iteration_history_tr_list[i], 'Total'))
    tgn.append(extract_history(iteration_history_gn_list[i], 'Total'))
    
c_TR_Xu_BPSK = np.array(ctr) # simulation 5
c_GN_Xu_BPSK = np.array(cgn) # simulation 3
Total_TR_Xu_BPSK = np.array(ttr) # simulation 5
Total_GN_Xu_BPSK = np.array(tgn) # simulation 3

#%% plot
idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]
idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
plt.figure(1)
# spectrum usage in simulation 3 and 5
plt.plot(np.arange(5, 55, 5), c_GN_Xu_BPSK.mean(axis=0)[idx], label='GN')
plt.plot(np.arange(5, 55, 5), c_TR_Xu_BPSK.mean(axis=0)[idx], label='TR')
plt.title('Spectrum Xu BPSK')
plt.legend()
plt.show()

plt.figure(2)
# number of regenerators in simulation 3 and 5
plt.plot(np.arange(5, 55, 5), Total_GN_Xu_BPSK.mean(axis=0)[idx], label='GN')
plt.plot(np.arange(5, 55, 5), Total_TR_Xu_BPSK.mean(axis=0)[idx], label='TR')
plt.title('#Regenerators Xu BPSK')
plt.legend()
plt.show()