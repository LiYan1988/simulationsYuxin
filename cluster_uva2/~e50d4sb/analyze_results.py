# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 16:25:16 2017

@author: yx4vf

analyze results from Rivanna
"""


from milp4 import *
import matplotlib.pyplot as plt

#%% running time
#total_time = []
#for i in range(4):
#    sn, iteration_history = read_data('e50d4sb_{}.pkl'.format(i))
#    print('Total time for simulation {}: {}'.format(i, sn.total_runtime/3600))
#    total_time.append(sn.total_runtime/3600)
#total_time = np.array(total_time)
#print('Mean running time: {}'.format(total_time.mean()))
#print('Max - Min running time: {}'.format(total_time.max()-total_time.min()))

#%%
colors = np.array([[0, 0.45, 0.74], [0.85, 0.33, 0.1], [0.93, 0.69, 0.13],
                   [0.49, 0.18, 0.50]])
hf1 = []
hf2 = []
hf3 = []

ctr = {}
cgn = {}
ttr = {}
tgn = {}
nnntr = {}
nnngn = {}
demands_solved = {}
stage_end = {}

demands_file = 'e50d4sb_0.csv'
demands = pd.read_csv(demands_file)
n_demands = demands.shape[0]
n_stages = int(np.ceil(n_demands/n_demands_per_stage))
n_demands_initial = n_demands-(n_stages-1)*n_demands_per_stage
n_demands_in_stage = [n_demands_initial+idx_stage*n_demands_per_stage
                      for idx_stage in range(n_stages)]

for i in range(4):
    sn, iteration_history = read_data('e50d4sb_{}.pkl'.format(i))
    demands_solved[i] = np.array([len(iteration_history[0][i]['demands_solved'])
        for i in iteration_history[0].keys()])
    stage_end[i] = [np.max(np.where(demands_solved[i]==j)[0])
        for j in n_demands_in_stage]

    # spectrum
    ctr[i] = extract_history(iteration_history[0], 'c')
    cgn[i] = extract_history(iteration_history[1], 'c')
    ctr[i] = [ctr[i][j] for j in stage_end[i]]
    cgn[i] = [cgn[i][j] for j in stage_end[i]]

    # total regeneraters
    ttr[i] = extract_history(iteration_history[0], 'Total')
    tgn[i] = extract_history(iteration_history[1], 'Total')
    ttr[i] = [ttr[i][j] for j in stage_end[i]]
    tgn[i] = [tgn[i][j] for j in stage_end[i]]

    # total circuit
    nnntr[i] = extract_history(iteration_history[0], 'NNN')
    nnngn[i] = extract_history(iteration_history[1], 'NNN')
    nnntr[i] = [sum(i.values()) for i in nnntr[i]]
    nnngn[i] = [sum(i.values()) for i in nnngn[i]]
    nnntr[i] = [nnntr[i][j] for j in stage_end[i]]
    nnngn[i] = [nnngn[i][j] for j in stage_end[i]]

ctrm = np.array([ctr[i] for i in ctr.keys()]).mean(axis=0)
cgnm = np.array([cgn[i] for i in cgn.keys()]).mean(axis=0)
plt.figure(1)
h1, = plt.plot(ctrm, label='TR spectrum', color=colors[0], linestyle='-.')
h2, = plt.plot(cgnm, label='GN spectrum', color=colors[1], linestyle='--')
plt.legend(handles=[h1, h2])

ttrm = np.array([ttr[i] for i in ttr.keys()]).mean(axis=0)
tgnm = np.array([tgn[i] for i in tgn.keys()]).mean(axis=0)
plt.figure(2)
h1, = plt.plot(ttrm, label='TR Total', color=colors[0], linestyle='-.')
h2, = plt.plot(tgnm, label='GN Total', color=colors[1], linestyle='--')
plt.legend(handles=[h1, h2])

nnntrm = np.array([nnntr[i] for i in nnntr.keys()]).mean(axis=0)
nnngnm = np.array([nnngn[i] for i in nnngn.keys()]).mean(axis=0)
plt.figure(3)
h1, = plt.plot(nnntrm, label='TR NNN', color=colors[0], linestyle='-.')
h2, = plt.plot(nnngnm, label='GN NNN', color=colors[1], linestyle='--')
plt.legend(handles=[h1, h2])

