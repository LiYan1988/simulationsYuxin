# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 18:44:46 2017

@author: yx4vf

test one sample
"""

from milp4_tmp import *
np.random.seed(0)

for k in range(20):
    network_cost = np.array([[INF, 5, INF, INF, INF, 6],
                             [5, INF, 4, INF, 5, 3],
                             [INF, 4, INF, 5, 3, INF],
                             [INF, INF, 5, INF, 6, INF],
                             [INF, 5, 3, 6, INF, 4],
                             [6, 3, INF, INF, 4, INF]])
    sn = Network(network_cost)
    n_demands = 20
    demands = sn.create_demands(n_demands, modulation='bpsk', low=40, high=100)
    demands.to_csv('6node-opt-demand{}.csv'.format(k))

    iteration_history_tr, iteration_history_gn = \
        sn.iterate(demands, miphint=True, mipfocus=1, method=-1, mipgap=0.001,
                   presolve=2, LogFile='log-opt-{}.log'.format(k),
                   max_added=half_n)

    iteration_history = (iteration_history_tr, iteration_history_gn)
    save_data('test_6node_opt{}.pkl'.format(k), iteration_history)

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

nnntr = extract_history(iteration_history_tr, 'NNN')
nnngn = extract_history(iteration_history_gn, 'NNN')
nnntr = [sum(i.values()) for i in nnntr]
nnngn = [sum(i.values()) for i in nnngn]

plt.figure(3)
plt.plot(nnntr)
plt.plot(nnngn)