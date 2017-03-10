# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 18:44:46 2017

@author: yx4vf

test one sample
"""

from milp4 import *

iteration_history = read_data('test2_milp4_cluster.pkl')


#%% 
ctr = extract_history(iteration_history[0], 'c')
ttr = extract_history(iteration_history[0], 'Total')
cgn = extract_history(iteration_history[1], 'c')
tgn = extract_history(iteration_history[1], 'Total')

import matplotlib.pyplot as plt
plt.figure(1)
plt.plot(ctr, label='TR')
plt.plot(cgn, label='GN')
plt.figure(2)
plt.plot(ttr, label='TR')
plt.plot(tgn, label='GN')

nnntr = extract_history(iteration_history[0], 'NNN')
nnngn = extract_history(iteration_history[1], 'NNN')
nnntr = [sum(i.values()) for i in nnntr]
nnngn = [sum(i.values()) for i in nnngn]

plt.figure(3)
plt.plot(nnntr)
plt.plot(nnngn)