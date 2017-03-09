# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 18:44:46 2017

@author: yx4vf

test one sample
"""

from milp2 import *
np.random.seed(0)

#network_cost = pd.read_csv('networkDT.csv', header=None)/100
network_cost = pd.read_csv('nsf-24nodes.csv', header=None, index_col=None)
network_cost = network_cost.as_matrix()
sn = Network(network_cost)
n_demands = 25
demands = sn.create_demands(n_demands, modulation='bpsk', low=40, high=100)
demands.to_csv('nsf24-demand1.csv')
#demands = demands.iloc[:5]

iteration_history_tr, iteration_history_gn = \
    sn.iterate(demands, random_state=0, miphint=True, mipfocus=1, 
               method=-1, mipgap=0.001, presolve=2)

# gurobi model cannot be save by pickle
models_gn = {}
models_tr = {}
for i in iteration_history_gn.keys():
    models_gn[i] = iteration_history_gn[i].pop('model', None)
    models_tr[i] = iteration_history_tr[i].pop('model', None)

iteration_history = (iteration_history_tr, iteration_history_gn)
save_data('nsf24-test1.pkl', iteration_history)

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

#objgn = [models_gn[i].ObjVal for i in models_gn.keys()]
#objtr = [models_tr[i].ObjVal for i in models_tr.keys()]
#
#plt.figure(3)
#plt.plot(objtr)
#plt.plot(objgn)

nnntr = extract_history(iteration_history_tr, 'NNN')
nnngn = extract_history(iteration_history_gn, 'NNN')
nnntr = [sum(i.values()) for i in nnntr]
nnngn = [sum(i.values()) for i in nnngn]

plt.figure(3)
plt.plot(nnntr)
plt.plot(nnngn)