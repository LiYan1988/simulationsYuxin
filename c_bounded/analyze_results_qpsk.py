# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 11:05:25 2017

@author: misWin

analyze simulation46xu data
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
    file_name = 'results/results46xu/simulation46xu_{}.pkl'.format(i)
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
    
c_TR_Xu_BPSK = np.array(ctr) # simulation 6
c_GN_Xu_BPSK = np.array(cgn) # simulation 4
Total_TR_Xu_BPSK = np.array(ttr) # simulation 6
Total_GN_Xu_BPSK = np.array(tgn) # simulation 4

#%% plot
idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]
idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
#idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
xaxis = np.arange(5, 55, 5)
plt.figure(1)
# spectrum usage in simulation 4 and 6
plt.plot(xaxis, c_GN_Xu_BPSK.mean(axis=0)[idx], label='GN')
plt.plot(xaxis, c_TR_Xu_BPSK.mean(axis=0)[idx], label='TR')
plt.title('Spectrum Xu QPSK')
plt.legend()
plt.show()

plt.figure(2)
# number of regenerators in simulation 4 and 6
plt.plot(xaxis, Total_GN_Xu_BPSK.mean(axis=0)[idx], label='GN')
plt.plot(xaxis, Total_TR_Xu_BPSK.mean(axis=0)[idx], label='TR')
plt.title('#Regenerators Xu QPSK')
plt.legend()
plt.show()