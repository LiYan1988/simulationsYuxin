# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 22:19:29 2017

@author: misWin
This is a template of python files for simulation 4 and 6 with Xu's algorithm on hebbe
"""


from milp2_xu import *
np.random.seed(0)

batch_id = 22
network_cost = pd.read_csv('nsf-24nodes.csv', header=None, index_col=None)
network_cost = network_cost.as_matrix()
sn = Network(network_cost, modulation='qpsk')
demands_file = 'simulation46xu_22.csv'
demands = pd.read_csv(demands_file)

iteration_history_tr, iteration_history_gn = \
    sn.iterate(demands, random_state=0, mipstart=True, mipfocus=1, 
               method=-1, mipgap=0.01)

# gurobi model instances cannot be save by pickle
models_gn = {}
models_tr = {}
for i in iteration_history_gn.keys():
    models_gn[i] = iteration_history_gn[i].pop('model', None)
    models_tr[i] = iteration_history_tr[i].pop('model', None)
    
iteration_history = (iteration_history_tr, iteration_history_gn)
output_file = 'simulation46xu_22.pkl'
save_data(output_file, iteration_history)