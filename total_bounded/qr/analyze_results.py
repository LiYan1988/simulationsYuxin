# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 11:08:44 2017

@author: yx4vf
"""

from milp4_1 import *
import matplotlib.pyplot as plt


tr_total = []
gn_total = []
tr_c = []
gn_c = []

demands_solved = {}
stage_end = {}
n_demands = 50
n_stages = int(np.ceil(n_demands/n_demands_per_stage))
n_demands_initial = n_demands-(n_stages-1)*n_demands_per_stage
n_demands_in_stage = [n_demands_initial+idx_stage*n_demands_per_stage
                      for idx_stage in range(n_stages)]

for i in range(5):
    sn, iteration_history = read_data('qr_{}.pkl'.format(i))
    demands_solved[i] = np.array([len(iteration_history[0][i]['demands_solved'])
        for i in iteration_history[0].keys()])
    stage_end[i] = [np.max(np.where(demands_solved[i]==j)[0])
        for j in n_demands_in_stage]

    iteration_history_tr = iteration_history[0]
    iteration_history_gn = iteration_history[1]
    tr_total.append(extract_history(iteration_history_tr, 'Total'))
    gn_total.append(extract_history(iteration_history_gn, 'Total'))
    tr_c.append(extract_history(iteration_history_tr, 'c'))
    gn_c.append(extract_history(iteration_history_gn, 'c'))

tr_total = np.array([[tr_total[i][j] for j in stage_end[i]] for i in range(5)])
gn_total = np.array([[gn_total[i][j] for j in stage_end[i]] for i in range(5)])
tr_c = np.array([[tr_c[i][j] for j in stage_end[i]] for i in range(5)])
gn_c = np.array([[gn_c[i][j] for j in stage_end[i]] for i in range(5)])

plt.figure(1)
plt.plot(tr_total.mean(axis=0), label='tr', linestyle='--')
plt.plot(gn_total.mean(axis=0), label='GN ', linestyle=':')
plt.legend()
plt.show()

plt.figure(2)
plt.plot(tr_c.mean(axis=0), label='TR proposed', linestyle='--')
plt.plot(gn_c.mean(axis=0), label='TR optimal', linestyle=':')
plt.legend()
plt.show()