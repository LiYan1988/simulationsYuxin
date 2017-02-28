# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 11:49:07 2017

@author: yx4vf

analyze simulation1 data
QPSK, TR and GN models, proposed algorithm
"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import gc

from milp2 import *

#%% read data
ctr = []
cgn = []
ttr = []
tgn = []
for i in range(30):
    print('Reading results of simulation {}'.format(i))
    try:
        file_name = 'results/results2_partial/simulation1_{}.pkl'.format(i)
        tr, gn = read_data(file_name)
        ctr.append(extract_history(tr, 'c'))
        cgn.append(extract_history(gn, 'c'))
        ttr.append(extract_history(tr, 'Total'))
        tgn.append(extract_history(gn, 'Total'))
        del tr, gn
        gc.collect()
    except:
        pass
    
#%% extract c and total    
c_TR_BPSK = np.array(ctr) # simulation 6
c_GN_BPSK = np.array(cgn) # simulation 4
Total_TR_BPSK = np.array(ttr) # simulation 6
Total_GN_BPSK = np.array(tgn) # simulation 4

#%% plot iteration history before all the simulations are done
plt.figure(1)
plt.plot(c_TR_BPSK.mean(axis=0), label='TR')
plt.plot(c_GN_BPSK.mean(axis=0), label='GN')
plt.legend()
plt.title('Spectrum')
plt.show()

plt.figure(2)
plt.plot(Total_TR_BPSK.mean(axis=0), label='TR')
plt.plot(Total_GN_BPSK.mean(axis=0), label='GN')
plt.legend()
plt.title('#Regenerators')
plt.show()

#%% plot
xaxis = np.arange(10, 110, 10)
idx = list(xaxis)
plt.figure(3)
# spectrum usage in simulation 4 and 6
plt.plot(xaxis, c_GN_BPSK.mean(axis=0)[idx], label='GN')
plt.plot(xaxis, c_TR_BPSK.mean(axis=0)[idx], label='TR')
plt.title('Spectrum SA BPSK')
plt.legend()
plt.show()

plt.figure(4)
# number of regenerators in simulation 4 and 6
plt.plot(xaxis, Total_GN_BPSK.mean(axis=0)[idx], label='GN')
plt.plot(xaxis, Total_TR_BPSK.mean(axis=0)[idx], label='TR')
plt.title('#Regenerators SA BPSK')
plt.legend()
plt.show()