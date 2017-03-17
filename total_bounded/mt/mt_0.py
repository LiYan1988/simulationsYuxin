# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 22:06:02 2017

@author: misWin

This is a template of python files for simulation algorithem on Rivanna
"""

from milp4_1 import *
seed = 224190
np.random.seed(seed)
kwargs = {'mipfocus': 1, 'miphint': True, 'presolve': 2, 'logfile': 'mt_0.log', 'mipgap': 0.001, 'threads': 4}

batch_id = 0
network_cost = pd.read_csv('nsf-24nodes.csv', header=None, index_col=None)
network_cost = network_cost.as_matrix()
sn = Network(network_cost, modulation='qpsk')
demands_file = 'mt_0.csv'
demands = pd.read_csv(demands_file)

iteration_history_tr, iteration_history_gn = sn.iterate(demands, **kwargs)

iteration_history = (iteration_history_tr, iteration_history_gn)
output_file = 'mt_0.pkl'
save_data(output_file, (sn, iteration_history))