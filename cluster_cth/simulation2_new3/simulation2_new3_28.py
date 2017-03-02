# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 20:10:26 2017

@author: misWin

This is a template of python files for simulation 2 on hebbe
"""

from milp2 import *
np.random.seed(0)

batch_id = 28
network_cost = pd.read_csv('nsf-24nodes.csv', header=None, index_col=None)
network_cost = network_cost.as_matrix()
sn = Network(network_cost, modulation='qpsk')
demands_file = 'simulation2_new3_28.csv'
demands = pd.read_csv(demands_file)

iteration_history_tr, iteration_history_gn = \
    sn.iterate(demands, random_state=0, mipstart=True, mipfocus=1, 
               method=-1, mipgap=0.001, Threads=20)

# gurobi model instances cannot be save by pickle
#models_gn = {}
#models_tr = {}
#for i in iteration_history_gn.keys():
#    models_gn[i] = iteration_history_gn[i].pop('model', None)
#    models_tr[i] = iteration_history_tr[i].pop('model', None)
    
iteration_history = (iteration_history_tr, iteration_history_gn)
output_file = 'simulation2_new3_28.pkl'
save_data(output_file, iteration_history)