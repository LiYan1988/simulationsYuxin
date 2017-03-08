# -*- coding: utf-8 -*-
"""
Created on Sat Feb 25 22:06:02 2017

@author: misWin

This is a template of python files for simulation 3 and 5 with Xu's algorithem on hebbe
"""

from milp3 import *
np.random.seed(0)

hostname = socket.gethostname()
memory = psutil.virtual_memory()
cpus = psutil.cpu_count()
node_info = (hostname, memory, cpus)
print('hostname: {}'.format(hostname))
print('memory: {}'.format(memory))
print('cpus: {}'.format(cpus))

batch_id = 0
network_cost = pd.read_csv('dt-14nodes.csv', header=None, index_col=None)
network_cost = network_cost.as_matrix()
sn = Network(network_cost, modulation='bpsk')
demands_file = '../demands/demands_14nodes_matlab_'+str(batch_id)+'.csv'
demands = pd.read_csv(demands_file).iloc[:10]

iteration_history = \
    sn.iterate(demands, random_state=0, mipstart=True, mipfocus=1, 
               method=-1, mipgap=0.001)    

output_file = 'output-GN-vs-TR-bpsk-nsf24'+str(batch_id)+'.pkl'
save_data(output_file, (iteration_history, node_info))
