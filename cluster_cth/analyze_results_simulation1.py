# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 11:50:14 2017

@author: misWin

analyze simulation1 data
BPSK, TR and GN models, proposed algorithm
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from milp2 import *

#%% read data
iteration_history_tr_list = []
iteration_history_gn_list = []
for i in range(30):
    print('Reading results of simulation {}'.format(i))
    try:
        file_name = 'results/results1_partial/simulation1_{}.pkl'.format(i)
        tmp1, tmp2 = read_data(file_name)
        iteration_history_tr_list.append(tmp1)
        iteration_history_gn_list.append(tmp2)
    except:
        pass
    
#%% extract c and total
ctr = []
cgn = []
ttr = []
tgn = []
for i in range(30):
    print('Extracting data from simulation {}'.format(i))
    try:
        ctr.append(extract_history(iteration_history_tr_list[i], 'c'))
        cgn.append(extract_history(iteration_history_gn_list[i], 'c'))
        ttr.append(extract_history(iteration_history_tr_list[i], 'Total'))
        tgn.append(extract_history(iteration_history_gn_list[i], 'Total'))
    except:
        pass
    
c_TR_Xu_BPSK = np.array(ctr) # simulation 6
c_GN_Xu_BPSK = np.array(cgn) # simulation 4
Total_TR_Xu_BPSK = np.array(ttr) # simulation 6
Total_GN_Xu_BPSK = np.array(tgn) # simulation 4

#%% plot iteration history before all the simulations are done
plt.figure(1)
plt.plot(c_TR_Xu_BPSK.mean(axis=0), label='TR')
plt.plot(c_GN_Xu_BPSK.mean(axis=0), label='GN')
plt.legend()
plt.title('Spectrum')
plt.show()

plt.figure(2)
plt.plot(Total_TR_Xu_BPSK.mean(axis=0), label='TR')
plt.plot(Total_GN_Xu_BPSK.mean(axis=0), label='GN')
plt.legend()
plt.title('#Regenerators')
plt.show()

#%% plot
#idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10]
#idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
##idx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
#xaxis = np.arange(5, 55, 5)
#plt.figure(1)
## spectrum usage in simulation 4 and 6
#plt.plot(xaxis, c_GN_Xu_BPSK.mean(axis=0)[idx], label='GN')
#plt.plot(xaxis, c_TR_Xu_BPSK.mean(axis=0)[idx], label='TR')
#plt.title('Spectrum Xu QPSK')
#plt.legend()
#plt.show()
#
#plt.figure(2)
## number of regenerators in simulation 4 and 6
#plt.plot(xaxis, Total_GN_Xu_BPSK.mean(axis=0)[idx], label='GN')
#plt.plot(xaxis, Total_TR_Xu_BPSK.mean(axis=0)[idx], label='TR')
#plt.title('#Regenerators Xu QPSK')
#plt.legend()
#plt.show()