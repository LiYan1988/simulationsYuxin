# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 22:06:02 2017

@author: misWin

This is a template of python files for simulation algorithem on Rivanna
"""

from milp4 import *
seed = 71860
np.random.seed(seed)
kwargs = {'logfile': 'bpsk_rest_7.log', 'mipgap': 0.001, 'miphint': True, 'threads': 4, 'mipfocus': 1, 'presolve': 2}

batch_id = 7
network_cost = pd.read_csv('nsf-24nodes.csv', header=None, index_col=None)
network_cost = network_cost.as_matrix()
sn = Network(network_cost, modulation='bpsk')
demands_file = 'bpsk_rest_7.csv'
demands = pd.read_csv(demands_file)

iteration_history_tr, iteration_history_gn = sn.iterate(demands, **kwargs)

iteration_history = (iteration_history_tr, iteration_history_gn)
output_file = 'bpsk_rest_7.pkl'
save_data(output_file, (sn, iteration_history))