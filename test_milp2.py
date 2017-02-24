# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 18:44:46 2017

@author: yx4vf
"""

from milp2 import *

network_cost = pd.read_csv('networkDT.csv', header=None)/100
network_cost = network_cost.as_matrix()
sn = Network(network_cost)
n_demands = 50
demands = sn.create_demands(n_demands, modulation='bpsk', low=40, high=100)
demands = demands.iloc[:30]

iteration_history_tr, iteration_history_gn = \
    sn.iterate(demands, random_state=0, mipstart=True, mipfocus=1, 
               timelimit=10, method=-1, mipgap=0.001, outputflag=1, 
               FeasibilityTol=1e-7, IntFeasTol=1e-7, OptimalityTol=1e-7)