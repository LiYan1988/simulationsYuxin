# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 22:06:02 2017

@author: misWin

This is a template of python files for simulation algorithem on Rivanna
"""

from milp4_1 import *
seed = 0
np.random.seed(seed)
kwargs = {'miphint':True, 'mipfocus':1, 'mipgap':0.001, 'presolve':2, 'logfile':'log_qpsk.log'}

batch_id = 0
network_cost = pd.read_csv('nsf-24nodes.csv', header=None, index_col=None)
network_cost = network_cost.as_matrix()
sn = Network(network_cost, modulation='qpsk')
demands_file = 'demands_template_'+str(batch_id)+'.csv'
demands = pd.read_csv(demands_file)

iteration_history_tr, iteration_history_gn = sn.iterate(demands, **kwargs)

iteration_history = (iteration_history_tr, iteration_history_gn)
output_file = 'output-GN-vs-TR-qpsk-nsf24'+str(batch_id)+'.pkl'
save_data(output_file, (sn, iteration_history))