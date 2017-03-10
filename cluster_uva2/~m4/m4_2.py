# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 22:06:02 2017

@author: misWin

This is a template of python files for simulation 3 and 5 with Xu's algorithem on hebbe
"""

from milp4 import *
seed = 398701
np.random.seed(seed)
kwargs = {'miphint': True, 'threads': 4, 'presolve': 2, 'mipgap': 0.001, 'mipfocus': 1, 'logfile': 'm4_2.log'}

batch_id = 2
network_cost = pd.read_csv('nsf-24nodes.csv', header=None, index_col=None)
network_cost = network_cost.as_matrix()
sn = Network(network_cost, modulation='bpsk')
demands_file = 'm4_2.csv'
demands = pd.read_csv(demands_file)

iteration_history_tr, iteration_history_gn = sn.iterate(demands, **kwargs)

iteration_history = (iteration_history_tr, iteration_history_gn)
output_file = 'm4_2.pkl'
save_data(output_file, (sn, iteration_history))